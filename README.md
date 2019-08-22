[![Circle CI](https://circleci.com/gh/onecodex/needletail.svg?style=shield&circle-token=65c2b7d87452dba5e8e3e967133311af478632a4)](https://circleci.com/gh/onecodex/needletail)
[![crates.io](https://img.shields.io/crates/v/needletail.svg)](https://crates.io/crates/needletail)

# Needletail

Needletail is a MIT-licensed, minimal-copying FASTA/FASTQ parser and _k_-mer processing library for Rust.

The goal is to write a fast *and* well-tested set of functions that more specialized bioinformatics programs can use.
Needletail's goal is to be as fast as the [readfq](https://github.com/lh3/readfq) C library at parsing FASTX files and much (i.e. 25 times) faster than equivalent Python implementations at _k_-mer counting.

# Example

```rust
extern crate needletail;
use std::env;
use std::fs::File;
use needletail::parse_sequences;

fn main() {
  let filename: String = env::args().nth(1).unwrap();

  let mut n_bases = 0;
  let mut n_valid_kmers = 0;
  parse_sequences(File::open(filename).expect("missing file"), |_| {}, |seq| {
    // seq.id is the name of the record
    // seq.seq is the base sequence
    // seq.qual is an optional quality score

    // keep track of the total number of bases
    n_bases += seq.seq.len();
    
    // keep track of the number of AAAA (or TTTT via canonicalization) in the 
    // file (normalize makes sure ever base is capitalized for comparison)
    for (_, kmer, _) in seq.normalize(false).kmers(4, true) {
      if kmer == b"AAAA" {
        n_valid_kmers += 1;
      }
    }
  }).expect("parsing failed");
  println!("There are {} bases in your file.", n_bases);
  println!("There are {} AAAAs in your file.", n_valid_kmers);
}
```

# Installation

Needletail requires `rust` and `cargo` to be installed.
Please use either your local package manager (`homebrew`, `apt-get`, `pacman`, etc) or install these via [rustup](https://www.rustup.rs/).

Once you have Rust set up, you can include needletail in your `Cargo.toml` file like:
```shell
[dependencies]
needletail = "^0.3.0"
```

To install needletail itself for development:
```shell
git clone https://github.com/onecodex/needletail
cargo test  # to run tests
```

# Getting Help

Questions are best directed as GitHub issues. We plan to add more documentation soon, but in the meantime "doc" comments are included in the source.

# Contributing

Please do! We're happy to discuss possible additions and/or accept pull requests.
