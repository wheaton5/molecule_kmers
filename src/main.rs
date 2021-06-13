#[macro_use]
extern crate clap;
extern crate fnv;
extern crate hashbrown;
extern crate flate2; 
extern crate itertools;
extern crate needletail;
extern crate rayon;
extern crate byteorder;
use rayon::prelude::*;
use needletail::{parse_fastx_file, Sequence, FastxReader};
use flate2::read::GzDecoder;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Read;
use itertools::izip;
use byteorder::{WriteBytesExt, LittleEndian};

use std::fs::File;
use std::io::{BufWriter, Write};

use hashbrown::{HashMap, HashSet};
use std::str;


use clap::{App};

fn main() {
    let params = load_params();
    eprintln!("load kmers");
    let (kmers, kmer_type) = load_kmers(&params);
    eprintln!("txg");
    if params.txg_r1s.len() > 0 {
        eprintln!("{:?}\t{:?}\t{:?}\t{:?}", params.txg_r1s, params.txg_trim_r1s, params.txg_r2s, params.txg_trim_r2s);
        process_txg(&params, &kmers);
    }
    if params.hic_r1s.len() > 0 {
        eprintln!("hic");
        process_hic(&params, &kmer_type);
    }
    if params.long_reads.len() > 0 {
        eprintln!("longreads");
        process_longreads(&params, &kmers);
    }
    if params.fasta != None {
        eprintln!("fasta");
        process_fasta(&params, &kmers);
    }
}



struct Kmers {
    paired_hets: HashMap<Vec<u8>, i32>,
    unpaired_hets: HashMap<Vec<u8>, i32>,
    homozygous: HashMap<Vec<u8>, i32>,
}

#[derive(Debug)]
#[derive(PartialEq)]
enum KMER_TYPE {
    PAIRED_HET,
    UNPAIRED_HET,
    HOM,
}

fn load_kmers(params: &Params) -> (HashMap<Vec<u8>, i32>, HashMap<Vec<u8>, (i32, KMER_TYPE)>) {
    let mut kmers: HashMap<Vec<u8>, i32> = HashMap::default();
    let mut kmer_id: i32 = 1; // no 0 kmer id as we are using sign for which of the pair
    let mut kmer_type: HashMap<Vec<u8>, (i32, KMER_TYPE)> = HashMap::new();
    if let Some(het_kmers) = &params.paired_kmers {
        let reader = File::open(het_kmers).expect("cannot open barcode file");
        let mut reader = BufReader::new(reader);
        let mut buf1 = vec![];
        let mut buf2 = vec![];
        let mut buf3 = vec![];
        let mut buf4 = vec![];
        loop {
            buf1.clear(); buf2.clear(); buf3.clear(); buf4.clear();
            let bytes1 = reader.read_until(b'\t', &mut buf1).expect("cannot read file");
            if bytes1 == 0 { break; } 
            
            let bytes2 = reader.read_until(b'\t', &mut buf2).expect("cannot read file");
            if bytes2 == 0 { break; } 
            
            let bytes3 = reader.read_until(b'\t', &mut buf3).expect("cannot read file");
            if bytes3 == 0 { break; }
            
            let tag = str::from_utf8(&buf3).unwrap();
            let first = str::from_utf8(&buf1).unwrap();
            if tag == "HOM\t" {
                kmer_type.insert(buf1[0..(bytes1-1)].to_vec(), (kmer_id, KMER_TYPE::HOM));
                //eprintln!("got HOM {}",first);
                kmers.insert(buf1[0..(bytes1-1)].to_vec(), kmer_id);
            } else if tag == "HET\t" {
                kmer_type.insert(buf1[0..(bytes1-1)].to_vec(), (kmer_id, KMER_TYPE::UNPAIRED_HET));
                //eprintln!("got HET {}",first);
                kmers.insert(buf1[0..(bytes1-1)].to_vec(), kmer_id);
            } else {
                //eprintln!("ok just paired het {}",tag);
                kmer_type.insert(buf1[0..(bytes1-1)].to_vec(), (kmer_id, KMER_TYPE::PAIRED_HET));
                kmer_type.insert(buf3[0..(bytes1-1)].to_vec(), (kmer_id+1, KMER_TYPE::PAIRED_HET));
                kmers.insert(buf1[0..(bytes1-1)].to_vec(), kmer_id);
                kmers.insert(buf3[0..(bytes3-1)].to_vec(), kmer_id+1);  
            }         
            let bytes4 = reader.read_until(b'\n', &mut buf4).expect("cannot read file");
            if bytes4 == 0 { break; } 
            
            kmer_id += 2;
        }
    }
    println!("{} kmers",kmers.len());
    (kmers, kmer_type)
}

fn process_fasta(params: &Params, kmer_ids: &HashMap<Vec<u8>, i32>) {
    let mut fasta = "".to_string();
    if let Some(f) = &params.fasta {
        fasta = f.to_string();
    }
    let mut reader = parse_fastx_file(&fasta).expect("invalid path/file");
    let mut current_readname = String::new();
    let mut linedex = 0;
    let writer = File::create(format!("{}/fasta_kmers.bin",params.output))
        .expect("Unable to create file");
    let mut writer = BufWriter::new(writer);
    while let Some(record) = reader.next() {
        let mut variants: Vec<i32> = Vec::new();
        let record = record.expect("invalid record");
        let seq = record.normalize(false);
        let rc = seq.reverse_complement();
        for (position, kmer, canonical) in seq.canonical_kmers(params.kmer_size, &rc) {
            if let Some(kmer_id) = kmer_ids.get(kmer) {
                 if canonical {
                    variants.push(-*kmer_id);
                 } else { variants.push(*kmer_id); }
                 variants.push(position as i32);
            }
        }
        handle_read(&variants, &mut writer);
    }
}

fn process_longreads(params: &Params, kmer_ids: &HashMap<Vec<u8>, i32>) {
    let reads = &params.long_reads;
    
    let mut to_iterate: Vec<(usize, String)> = Vec::new();
    for (filenum, read) in reads.iter().enumerate() {
        to_iterate.push((filenum, read.to_string()));
    }
    let f = File::create(format!("{}/ccs.fofn", params.output)).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for (filenum, _r) in to_iterate.iter() {
        f.write_all(format!("{}/molecules_longreads_{:03}.bin\n", params.output, filenum).as_bytes()).expect("failed to write");
    }
    to_iterate.par_iter().for_each(|(filenum, read_file)|  {
        eprintln!("looking at file {}", read_file);
        let mut reader = parse_fastx_file(&read_file).expect("invalid path/file");
        let writer = File::create(format!("{}/molecules_longreads_{:03}.bin",params.output,filenum))
            .expect("Unable to create file");
        let mut writer = BufWriter::new(writer);
        while let Some(record) = reader.next() {
            let record = record.expect("invalid record");
            let seq = record.normalize(false);
            let rc = seq.reverse_complement();
            let mut variants: Vec<i32> = Vec::new();
            for (position, kmer, canonical) in seq.canonical_kmers(params.kmer_size, &rc) {
                if let Some(kmer_id) = kmer_ids.get(kmer) {
                    if canonical {
                        variants.push(-*kmer_id);
                    } else { variants.push(*kmer_id); }
                    variants.push(position as i32);
                }
            }
            handle_read(&variants, &mut writer);
        }
        
    });
}

fn handle_read(vars: &Vec<i32>, writer: &mut Write) {
    let mut result: Vec<u8> = Vec::new();
    for kmer in vars {
        result.write_i32::<LittleEndian>(*kmer).expect("write fail"); 
    } 
    result.write_i32::<LittleEndian>(0).expect("buffer fill fail");
    writer.write_all(&result).expect("write fail");
}

fn process_hic(params: &Params, kmer_type: &HashMap<Vec<u8>, (i32, KMER_TYPE)>) {
    let read1s = &params.hic_r1s;
    let read2s = &params.hic_r2s;
    let mut to_iterate: Vec<(usize, String, String)> = Vec::new();
    for (filenum, r1, r2) in izip!(0..read1s.len(), read1s, read2s) {
        to_iterate.push((filenum, r1.to_string(), r2.to_string()));
    }
    let f = File::create(format!("{}/hic.fofn", params.output)).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for (filenum, _r1, _r2) in to_iterate.iter() {
        f.write_all(format!("{}/molecules_hic_{:03}.bin\n", params.output, filenum).as_bytes()).expect("failed to write");
    }
    to_iterate.par_iter().for_each(|(filenum, r1_file, r2_file)| {
        let mut r1_reader = parse_fastx_file(&r1_file).expect("invalid path/file");
        let mut r2_reader = parse_fastx_file(&r2_file).expect("invalid path/file");

        let writer = File::create(format!("{}/molecules_hic_{:03}.bin",params.output,filenum))
            .expect("Unable to create file");
        let mut writer = BufWriter::new(writer);
        
        let mut vars1: Vec<i32> = Vec::new();
        let mut vars2: Vec<i32> = Vec::new();
        let mut result: Vec<u8> = Vec::new();
        while let Some(r1_record) = r1_reader.next() {
             if let Some(r2_record) = r2_reader.next() {
                let r1_rec = r1_record.expect("invalid record");
                let r1_seq = r1_rec.normalize(false);
                let r2_rec = r2_record.expect("invalid record");
                let r2_seq = r2_rec.normalize(false);
                let rc1 = r1_seq.reverse_complement();
                let rc2 = r2_seq.reverse_complement();
                for (_, kmer, _) in r1_seq.canonical_kmers(params.kmer_size, &rc1) {
                    if let Some((kmer_id, kmer_type)) = kmer_type.get(kmer) {
                        if *kmer_type == KMER_TYPE::PAIRED_HET {
                            vars1.push(*kmer_id);
                        } 
                    }
                }
                for (_, kmer, _) in r2_seq.canonical_kmers(params.kmer_size, &rc2) {
                    if let Some((kmer_id, kmer_type)) = kmer_type.get(kmer) {
                        if *kmer_type == KMER_TYPE::PAIRED_HET {
                            let mut already_has = false;
                            for var in vars1.iter() {
                                if var.abs() == kmer_id.abs() { already_has = true; }
                            }
                            if !already_has {
                                vars2.push(*kmer_id);
                            }
                        } 
                    }
                }
            }
            if vars1.len() > 0 && vars2.len() > 0 {
                for var in vars1.iter().chain(vars2.iter()) {
                    result.write_i32::<LittleEndian>(*var).expect("buffer fill fail");
                }
                result.write_i32::<LittleEndian>(0).expect("buffer fill fail");
                writer.write_all(&result).expect("write fail");
            }
            vars1.clear();
            vars2.clear();
            result.clear();
        }
    });
    //write_1_zero(); 
    //std::io::stdout().flush().expect("flush failed");
    //eprintln!("{} hic molecules written",molecules_written);
}

fn process_txg(params: &Params, kmer_ids: &HashMap<Vec<u8>, i32>)  {
    let barcodes = load_whitelist(params.txg_barcodes.as_ref().unwrap());
    let read1s = &params.txg_r1s;
    let read2s = &params.txg_r2s;
    let read1_trims = &params.txg_trim_r1s;
    let read2_trims = &params.txg_trim_r2s;
    let mut to_iterate: Vec<(usize, String, usize, String, usize)> = Vec::new();

    for (filenum, r1_file, r1_trim, r2_file, r2_trim) in 
        izip!(0..read1s.len(), read1s, read1_trims, read2s, read2_trims) {
        to_iterate.push((filenum, r1_file.to_string(), *r1_trim, r2_file.to_string(), *r2_trim));
    }
    let f = File::create(format!("{}/txg.fofn", params.output)).expect("Unable to create file");
    let mut f = BufWriter::new(f);
    for (filenum, r1, _, r2, _) in to_iterate.iter() {
        f.write_all(format!("{}/molecules_txg_{:03}.bin\n", params.output, filenum).as_bytes()).expect("failed to write");
    }
    
    to_iterate.par_iter().for_each(|(filenum, r1_file, r1_trim, r2_file, r2_trim)| {
        let mut r1_reader = parse_fastx_file(&r1_file).expect("invalid path/file");
        let mut r2_reader = parse_fastx_file(&r2_file).expect("invalid path/file");
        let mut result: Vec<u8> = Vec::new();
        let r1_trim = *r1_trim;
        let r2_trim = *r2_trim;
        let writer = File::create(format!("{}/molecules_txg_{:03}.bin",params.output,filenum))
            .expect("Unable to create file");
        let mut writer = BufWriter::new(writer);
        while let Some(r1_record) = r1_reader.next() {
            if let Some(r2_record) = r2_reader.next() {
                let r1_rec = r1_record.expect("invalid record");
                let r1_seq = r1_rec.normalize(false);
                if let Some(barcode_id) = barcodes.get(&r1_rec.sequence()[0..16]) {
                    let r2_rec = r2_record.expect("invalid record");
                    let r2_seq = r2_rec.normalize(false);
                    let rc1 = r1_seq.reverse_complement();
                    let rc2 = r2_seq.reverse_complement();
                    for (index, kmer, _) in r1_seq.canonical_kmers(params.kmer_size, &rc1) {
                        if index > r1_trim {
                            if let Some(kmer_id) = kmer_ids.get(kmer) {
                                result.write_i32::<LittleEndian>(*barcode_id).expect("buffer fill fail");
                                result.write_i32::<LittleEndian>(*kmer_id).expect("buffer fail");
                                writer.write_all(&result).expect("write fail");
                                result.clear();
                            }
                        }
                    }
                    for (index, kmer, _) in r2_seq.canonical_kmers(params.kmer_size, &rc2) {
                        if let Some(kmer_id) = kmer_ids.get(kmer) {
                            if index > r2_trim {
                                result.write_i32::<LittleEndian>(*barcode_id).expect("buffer fill fail");
                                result.write_i32::<LittleEndian>(*kmer_id).expect("buffer fail");
                                writer.write_all(&result).expect("write fail");
                                result.clear();
                            }
                        }
                    }
                }
            }
        }
    });
    //write_2_zeros(); // delimiter
    //std::io::stdout().flush().expect("flush failed");
}

fn load_whitelist(barcode_file: &String) -> HashMap<Vec<u8>, i32> {
    let mut buf: Vec<u8> = vec![];
    let mut barcodes: HashMap<Vec<u8>, i32> = HashMap::with_capacity(20000000);
    let mut barcode_index: i32 = 1;
    let reader = File::open(barcode_file).expect("cannot open barcode file");
    let mut reader = BufReader::new(reader);
    loop {
        let bytes = reader.read_until(b'\n', &mut buf).expect("cannot read file");
        if bytes == 0 { break; }
        barcodes.insert(buf[0..(bytes-1)].to_owned(), barcode_index);
        buf.clear();
        barcode_index += 1;
    }
    barcodes 
}

pub fn get_reader(filename: String) -> BufReader<Box<dyn Read>> {
    let filetype: Vec<&str> = filename.split(".").collect();
    let filetype = filetype[filetype.len()-1];
    let file = match File::open(filename.clone()) {
        Ok(file) => file,
        Err(error) => panic!("There was a problem opening the file: {:?}", error),
    };
    let reader: Box<dyn Read> = match filetype { 
        "gz" => Box::new(GzDecoder::new(file)), 
        _ => Box::new(file),
    }; 
    BufReader::new(reader)
}


#[derive(Clone)]
struct Params {
    paired_kmers: Option<String>,
    unpaired_kmers: Option<String>,
    txg_r1s: Vec<String>,
    txg_r2s: Vec<String>,
    txg_trim_r1s: Vec<usize>,
    txg_trim_r2s: Vec<usize>,
    txg_barcodes: Option<String>,
    long_reads: Vec<String>,
    hic_r1s: Vec<String>,
    hic_r2s: Vec<String>,
    output: String,
    kmer_size: u8,
    fasta: Option<String>,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let paired_kmers = match params.value_of("paired_kmers") {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    let unpaired_kmers = match params.value_of("unpaired_kmers") {
        Some(x) => Some(x.to_string()),
        None => None,
    };

    let mut txg_r1s: Vec<String> = Vec::new();
    let mut txg_r2s: Vec<String> = Vec::new();
    let mut txg_trim_r1s: Vec<usize> = Vec::new();
    let mut txg_trim_r2s: Vec<usize> = Vec::new();
    match params.value_of("txg_reads") {
        Some(txg_fofn) => {
            let f = File::open(txg_fofn).expect("Unable to open txg fofn");
            let f = BufReader::new(f);

            for (linedex, line) in f.lines().enumerate() {
                let line = line.expect("Unable to read txg fofn line");
                let vec: Vec<&str> = line.split_whitespace().collect();
                let trim = vec[1].to_string().parse::<usize>().unwrap();
                if linedex % 2 == 0 {
                    txg_r1s.push(vec[0].to_string());
                    txg_trim_r1s.push(trim);
                } else {
                    txg_r2s.push(vec[0].to_string());
                    txg_trim_r2s.push(trim);
                }
            }
        },
        None => (),
    }

    let mut hic_r1s: Vec<String> = Vec::new();
    let mut hic_r2s: Vec<String> = Vec::new();
    match params.value_of("hic_reads") {
        Some(hic_fofn) => {
            let f = File::open(hic_fofn).expect("Unable to open hic fofn");
            let f = BufReader::new(f);

            for (linedex, line) in f.lines().enumerate() {
                let line = line.expect("Unable to read hic fofn line");
                if linedex % 2 == 0 {
                    hic_r1s.push(line.to_string());
                } else {
                    hic_r2s.push(line.to_string());
                }
            }
        },
        None => (),
    }

    let txg_barcodes = match params.value_of("txg_barcodes") {
        Some(x) => Some(x.to_string()),
        None => {assert!(txg_r1s.len() == 0, "if you specify txg_files, you must also specify txg_barcodes");  None},
    };

    let mut long_reads: Vec<String> = Vec::new();
    match params.value_of("long_reads") {
        Some(ccs_fofn) => {
            let f = File::open(ccs_fofn).expect("Unable to open long read fofn");
            let f = BufReader::new(f);

            for line in f.lines() {
                let text = line.expect("Unable to read long read fofn line");
                let toks = text.split_whitespace();
                let toks: Vec<&str> = toks.collect();
                long_reads.push(toks[0].to_string());
            }
        },
        None => (),
    }

    let output = params.value_of("output").unwrap();

    let kmer_size = params.value_of("kmer_size").unwrap();
    let kmer_size: u8 = kmer_size.to_string().parse::<u8>().unwrap();

    let fasta: Option<String> = match params.value_of("fasta") {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    

    Params{
        fasta: fasta,
        output: output.to_string(),
        kmer_size: kmer_size,
        paired_kmers: paired_kmers,
        unpaired_kmers: unpaired_kmers,
        txg_r1s: txg_r1s,
        txg_r2s: txg_r2s,
        txg_trim_r1s: txg_trim_r1s, 
        txg_trim_r2s: txg_trim_r2s,
        txg_barcodes: txg_barcodes,
        long_reads: long_reads,
        hic_r1s: hic_r1s,
        hic_r2s: hic_r2s,
    }   
}

fn write_1_zero() {
    let mut result: Vec<u8> = Vec::new();
    result.write_i32::<LittleEndian>(0).expect("fail");
    std::io::stdout().lock().write_all(&result).expect("fail");
}

fn write_2_zeros() {
    let mut result: Vec<u8> = Vec::new();
    result.write_i32::<LittleEndian>(0).expect("fail");
    result.write_i32::<LittleEndian>(0).expect("fail");
    std::io::stdout().lock().write_all(&result).expect("fail");
}
