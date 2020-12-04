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
use needletail::Sequence;
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
//use fnv::FnvHasher;
//use std::hash::BuildHasherDefault;
//type FnvHashMap<T,V> = HashMap<T,V, BuildHasherDefault<FnvHasher>>;
//type FnvHashSet<T> = HashSet<T, BuildHasherDefault<FnvHasher>>;
//use hashbrown::hash_map::DefaultHashBuilder;

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

fn process_longreads(params: &Params, kmer_ids: &HashMap<Vec<u8>, i32>)  {
    let reads = &params.long_reads;
    
    //write_2_zeros(); // delimiter, this file format isnt great
    let mut to_iterate: Vec<(usize, String)> = Vec::new();
    for (filenum, read) in reads.iter().enumerate() {
        to_iterate.push((filenum, read.to_string()));
    }
    to_iterate.par_iter().for_each(|(filenum, read_file)|  {
        let mut reader = get_reader(read_file.to_string());
        let mut current_readname = String::new();
        let mut buf = vec![];
        let mut linedex = 0;
        //let mut molecule_vars: Vec<Vec<i32>> = Vec::new();

        let writer = File::create(format!("{}/molecules_longreads_{:03}.bin",params.output,filenum))
            .expect("Unable to create file");
        let mut writer = BufWriter::new(writer);
        //let mut vars: HashMap<i32, [u8; 2]> = HashMap::new(); // going to count how many of each allele across subreads to do consensus

        loop {
            buf.clear();
            let bytes = reader.read_until(b'\n', &mut buf).expect("cannot read longread fastq");
            if bytes == 0 {  break; } 
            if linedex % 4 == 1 {
                let mut variants: Vec<i32> = Vec::new();
                let seq = &buf.sequence(); 
                let seq = seq.normalize(false);
                let rc = seq.reverse_complement();
                for (position, kmer, canonical) in seq.canonical_kmers(params.kmer_size, &rc) {
                    if let Some(kmer_id) = kmer_ids.get(kmer) {
                        //let counts = vars.entry(kmer_id.abs()).or_insert([0u8; 2]);
                        //if kmer_id < &0 { counts[1] += 1; } else { counts[0] += 1; }
                        if canonical {
                            variants.push(-*kmer_id);
                        } else { variants.push(*kmer_id); }
                        variants.push(position as i32);
                    }
                }
                handle_subreads(&variants, &mut writer);
            }
            linedex += 1; buf.clear();
        }        
    });
}

//fn handle_subreads(vars: &HashMap<i32, [u8; 2]>) {
fn handle_subreads(vars: &Vec<i32>, writer: &mut Write) {
    let mut result: Vec<u8> = Vec::new();
    //let mut to_output: Vec<i32> = Vec::new();
    /*
    for (kmer_id, counts) in vars.iter() {
        if counts[0] > 3*counts[1] {
            to_output.push(*kmer_id);
        } else if counts[1] > 3*counts[0] {
            to_output.push(-kmer_id);
        }
    }
    */
    //for kmer_id in vars.iter() {
    //    to_output.push(*kmer_id);
    //}
    //if to_output.len() > 1 { // write 1 entry for every read (multiple subreads per read) even if there are no kmers... will make figuring out which 
        for kmer in vars {
           result.write_i32::<LittleEndian>(*kmer).expect("write fail"); 
        } 
        result.write_i32::<LittleEndian>(0).expect("buffer fill fail");
        //std::io::stdout().lock().write_all(&result).expect("fail");
        //write_1_zero();
        writer.write_all(&result).expect("write fail");
    //}
}

fn process_hic(params: &Params, kmer_type: &HashMap<Vec<u8>, (i32, KMER_TYPE)>) {
    let read1s = &params.hic_r1s;
    let read2s = &params.hic_r2s;
    let mut to_iterate: Vec<(usize, String, String)> = Vec::new();
    for (filenum, r1, r2) in izip!(0..read1s.len(), read1s, read2s) {
        to_iterate.push((filenum, r1.to_string(), r2.to_string()));
    }
    to_iterate.par_iter().for_each(|(filenum, r1_file, r2_file)| {
        let mut r1_reader = get_reader(r1_file.to_string());
        let mut r2_reader = get_reader(r2_file.to_string()); 
        let mut buf1 = vec![];
        let mut buf2 = vec![];

        let writer = File::create(format!("{}/molecules_hic_{:03}.bin",params.output,filenum))
            .expect("Unable to create file");
        let mut writer = BufWriter::new(writer);
        
        let mut vars: Vec<i32> = Vec::new();
        let mut result: Vec<u8> = Vec::new();
        let mut linedex = 0;
        loop {
            let bytes = r1_reader.read_until(b'\n', &mut buf1).expect("cannot read r1file");
            let bytes2 = r2_reader.read_until(b'\n', &mut buf2).expect("cannot read r2file");
            if bytes == 0 || bytes2 == 0 { break; }
            if linedex % 4 != 1 { 
                linedex += 1;
                buf1.clear();
                buf2.clear();
                continue; 
            }
            let r1_sequence = &buf1.sequence();
            let r1_sequence = r1_sequence.normalize(false);
            let r1_rc = r1_sequence.reverse_complement();
            for (_, kmer, _) in r1_sequence.canonical_kmers(params.kmer_size, &r1_rc) {
                if let Some((kmer_id, kmer_type)) = kmer_type.get(kmer) {
                    if *kmer_type == KMER_TYPE::PAIRED_HET {
                        vars.push(*kmer_id);
                    } 
                }
            }
            let r2_sequence = &buf2.sequence();
            let r2_sequence = r2_sequence.normalize(false);
            let r2_rc = r2_sequence.reverse_complement();
            for (_, kmer, _) in r2_sequence.canonical_kmers(params.kmer_size, &r2_rc) {
                if let Some((kmer_id, kmer_type)) = kmer_type.get(kmer) {
                    if *kmer_type == KMER_TYPE::PAIRED_HET {
                        vars.push(*kmer_id);
                    } 
                }
            }
            if vars.len() > 1 {
                for index1 in 0..vars.len() {
                    result.write_i32::<LittleEndian>(vars[index1]).expect("buffer fill fail");
                }
                result.write_i32::<LittleEndian>(0).expect("buffer fill fail");
                //std::io::stdout().lock().write_all(&result).expect("fail");
                writer.write_all(&result).expect("write fail");
            }
            
            result.clear();
            linedex += 1;
            vars.clear();
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
    to_iterate.par_iter().for_each(|(filenum, r1_file, r1_trim, r2_file, r2_trim)| {
        // get readers
        let mut r1_reader = get_reader(r1_file.to_string());
        let mut r2_reader = get_reader(r2_file.to_string());
        let mut buf1 = vec![];
        let mut buf2 = vec![];

        // get writer
        let writer = File::create(format!("{}/molecules_txg_{:03}.bin",params.output,filenum))
            .expect("Unable to create file");
        let mut writer = BufWriter::new(writer);
        let mut linedex = 0;
        let mut result: Vec<u8> = Vec::new();
        loop {
            let bytes = r1_reader.read_until(b'\n', &mut buf1).expect("cannot read r1file");
            let bytes2 = r2_reader.read_until(b'\n', &mut buf2).expect("cannot read r2file");
            if bytes == 0 || bytes2 == 0 { break; }
            if linedex % 4 != 1 { 
                linedex += 1;
                buf1.clear();
                buf2.clear();
                continue; 
            }
            if let Some(barcode_id) = barcodes.get(&buf1[0..16]) {
                let r1_sequence = &buf1[(16+r1_trim)..].sequence();
                let r1_sequence = r1_sequence.normalize(false);
                let r1_rc = r1_sequence.reverse_complement();
                for (_, kmer, _) in r1_sequence.canonical_kmers(params.kmer_size, &r1_rc) {
                //r1_sequence.canonical_kmers(21, &r1_rc).collect::<Vec<(usize, &[u8], bool)>>().into_par_iter().for_each(|(_, kmer, _)| {
                    if let Some(kmer_id) = kmer_ids.get(kmer) {
                        //eprintln!("txg kmer of type {:?}",kmer_type.get(kmer).unwrap());
                        //eprintln!("{}\t{}", barcode_id, std::str::from_utf8(&kmer).unwrap()); 
                        result.write_i32::<LittleEndian>(*barcode_id).expect("buffer fill fail");
                        result.write_i32::<LittleEndian>(*kmer_id).expect("buffer fail");
                        //std::io::stdout().lock().write_all(&result).expect("fail");
                        writer.write_all(&result).expect("write fail");
                        result.clear();
                    }
                }
                let r2_sequence = &buf2[*r2_trim..].sequence();
                let r2_sequence = r2_sequence.normalize(false);
                let r2_rc = r2_sequence.reverse_complement();
                for (_, kmer, _) in r2_sequence.canonical_kmers(params.kmer_size, &r2_rc) {
                //r2_sequence.canonical_kmers(21, &r2_rc).collect::<Vec<(usize, &[u8], bool)>>().into_par_iter().for_each(|(_, kmer, _)| {
                    if let Some(kmer_id) = kmer_ids.get(kmer) {
                        //println!("{}\t{}", barcode_id, kmer_id);//std::str::from_utf8(&kmer).unwrap()); 
                        let mut result: Vec<u8> = Vec::new();
                        result.write_i32::<LittleEndian>(*barcode_id).expect("buffer fill fail");
                        result.write_i32::<LittleEndian>(*kmer_id).expect("buffer fill fail");
                        //std::io::stdout().lock().write_all(&result).expect("fail");
                        writer.write_all(&result).expect("write fail");
                        result.clear();
                    }
                }
            }
            buf1.clear();
            buf2.clear();
            linedex += 1;
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

fn get_reader(filename: String) -> BufReader<Box<dyn Read>> {
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
    let txg_r1s_tmp = match params.values_of("txg_r1s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_r1s: Vec<String> = Vec::new();
    for x in txg_r1s_tmp { txg_r1s.push(x.to_string()); }
    let txg_r2s_tmp = match params.values_of("txg_r2s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_r2s: Vec<String> = Vec::new();
    for x in txg_r2s_tmp { txg_r2s.push(x.to_string()); }
    let txg_trim_tmp = match params.values_of("txg_trim_r1s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_trim_r1s: Vec<usize> = Vec::new();
    for trim in txg_trim_tmp {
        txg_trim_r1s.push(trim.parse::<usize>().unwrap());
    }
    let txg_trim_tmp = match params.values_of("txg_trim_r2s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut txg_trim_r2s: Vec<usize> = Vec::new();
    for trim in txg_trim_tmp {
        txg_trim_r2s.push(trim.parse::<usize>().unwrap());
    }
    assert!(txg_trim_r1s.len() == txg_r1s.len(), "If you specify txg_r1s, must also specify trim length for each file (txg_trim_r1s) {} {}",txg_r1s.len(), txg_trim_r1s.len());
    assert!(txg_trim_r2s.len() == txg_r2s.len(), "If you specify txg_r2s, must also specify trim length for each file (txg_trim_r2s) {} {}",txg_r2s.len(), txg_trim_r2s.len());
    let hic_r1s_tmp = match params.values_of("hic_r1s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut hic_r1s: Vec<String> = Vec::new();
    for x in hic_r1s_tmp { hic_r1s.push(x.to_string()); }
    let hic_r2s_tmp = match params.values_of("hic_r2s") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut hic_r2s: Vec<String> = Vec::new();
    for x in hic_r2s_tmp { hic_r2s.push(x.to_string()); }
    assert!(hic_r2s.len() == hic_r1s.len(), "If you specify hic_r1s or hic_r2s, must also specify the same number of the other.");

    let txg_barcodes = match params.value_of("txg_barcodes") {
        Some(x) => Some(x.to_string()),
        None => {assert!(txg_r1s.len() == 0, "if you specify txg_files, you must also specify txg_barcodes");  None},
    };
    let long_reads_tmp = match params.values_of("long_reads") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut long_reads: Vec<String> = Vec::new();
    for x in long_reads_tmp { long_reads.push(x.to_string()); }
    let output = params.value_of("output").unwrap();

    let kmer_size = params.value_of("kmer_size").unwrap();
    let kmer_size: u8 = kmer_size.to_string().parse::<u8>().unwrap();

    Params{
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
