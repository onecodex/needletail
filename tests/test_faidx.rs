use std::collections::HashMap;
use std::str;
use std::io::{BufReader, Seek, SeekFrom};
use needletail::FastxReader;
// use needletail::FastxReader;
use needletail::parser::FastaReader; // Only call FastaReader for `seek`
use needletail::parser::IndexedReader;


#[test]
fn test_faidx() {
    println!("test_faidx");
    let mut file_reader = IndexedReader::from_path("tests/data/test_index.fa").unwrap();
    file_reader.seek(0).unwrap(); // for record1
    let record =
        if let Some(record) = file_reader.reader.next() { record }
        else { panic!("invalid record!!!!") };
    let rec_seq = record.expect("invalid record");
    let id = rec_seq.id();
    let id_str = str::from_utf8(&*id).unwrap();
    let seq = &rec_seq.seq()[2..4];
    let seq_str = str::from_utf8(&*seq).unwrap();
    println!("Read Record1: {}:{}", id_str,seq_str);
    file_reader.seek(33).unwrap(); // for record2
    let record =
        if let Some(record) = file_reader.reader.next() { record }
        else { panic!("invalid record!!!!") };
    let rec_seq = record.expect("invalid record");
    let id = rec_seq.id();
    let id_str = str::from_utf8(&*id).unwrap();
    let seq = rec_seq.seq();
    let seq_str = str::from_utf8(&*seq).unwrap();
    println!("Read Record2: {}:{}", id_str,seq_str);
}