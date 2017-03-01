extern crate needletail;

use std::env;
use needletail::{fastx, kmer};


fn count_bases(filename: String) {
    let mut n_bases = 0;
    fastx::fastx_cli(&filename[..], |_| (), |seq| {
        n_bases += seq.1.len();
    }).unwrap();
    println!("{:?}", n_bases);
}


fn main() {
    let filename: String = env::args().nth(1).unwrap();

    count_bases(filename);
}
