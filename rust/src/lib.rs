#![crate_name = "needletail"]

use std::str;
use std::string;
use std::path::Path;
use std::fs::File;
use std::io::Read;
use std::borrow::Cow;

pub mod fastx;
pub mod kmer;

fn printcow(ckmer: Cow<[u8]>) {
    match ckmer {
        Cow::Borrowed(i) => {
            println!("ref!{}!", str::from_utf8(i).unwrap());
        },
        Cow::Owned(i) => {
            println!("own!{}!", String::from_utf8(i).unwrap());
        }
    }
}


fn main() {
//    let mut fr = fastx::FASTXReader::new(b"@ta\nA\n+t\nB");
//    println!("{}", fr.reset_eof());
//    let results = fr.next().unwrap();
//    println!("{}", results.0.unwrap());
//    printb(results.1);
//    printb(results.2.unwrap());
//    println!("{}", fr.reset_eof());

//    println!("Canonicalizing AAGT");
//    let mut cowstring = kmer::canonical(Cow::Borrowed(b"AAGT"));
//    printcow(cowstring);
//    println!("Canonicalizing TAGT");
//    cowstring = kmer::canonical(Cow::Borrowed(b"TAGT"));
//    printcow(cowstring);

    for kmer in kmer::Kmer::new(b"AGCA\nCT", 2, true) { //.filter(|kmer| !kmer.contains(&('N' as u8))) {
        println!("??");
        printcow(kmer.clone());
        printcow(kmer::canonical(kmer));
    }

//    let path = Path::new("./test.fa");
//    let mut file = match File::open(&path) {
//        Err(e) => panic!("Error opening file: {}", e),
//        Ok(f) => f
//    };
//
//    let mut buffer = Vec::new();
//    file.read_to_end(&mut buffer);
//    let contigs = fastx::FASTXReader::new(&buffer);
//    for ctg in contigs {
//        println!("{}", ctg.0.unwrap());
//        printb(ctg.1);
//    }
}
