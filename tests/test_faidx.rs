use std::collections::HashMap;
use std::str;
use std::io::{BufReader, Seek, SeekFrom};
use needletail::FastxReader;
// use needletail::FastxReader;
use needletail::parser::FastaReader; // Only call FastaReader for `seek`
use needletail::parser::IndexedReader;


#[test]
fn test_faidx() {
    println!("test faidx");
    let mut file_reader = IndexedReader::from_path("tests/data/test_index.fa").unwrap();
    let seq = file_reader.test_fetch("xxx", 2, 8);
    let seq_str = str::from_utf8(&*seq).unwrap();
    println!("Read Seq: \n{}", seq_str);

    // let rec_seq = file_reader.fetch(&"yyy".to_string()).unwrap(); // for record1
    // let id = rec_seq.id();
    // let id_str = str::from_utf8(&*id).unwrap();
    // let seq = &rec_seq.seq();
    // let seq_str = str::from_utf8(&*seq).unwrap();
    // println!("Read Record1: {}:{}", id_str,seq_str);
    // let rec_seq = file_reader.fetch(&"zzz".to_string()); // for record1
    // let id = rec_seq.id();
    // let id_str = str::from_utf8(&*id).unwrap();
    // let seq = rec_seq.seq();
    // let seq_str = str::from_utf8(&*seq).unwrap();
    // println!("Read Record2: {}:{}", id_str,seq_str);
}