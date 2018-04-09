[![Circle CI](https://circleci.com/gh/onecodex/needletail.svg?style=shield&circle-token=65c2b7d87452dba5e8e3e967133311af478632a4)](https://circleci.com/gh/onecodex/needletail)

# Needletail

Needletail is a MIT-licensed, minimal-copying FASTA/FASTQ parser and _k_-mer processing library for Rust.

The goal is to write a fast *and* well-tested set of functions that more specialized bioinformatics programs can use.
Needletail's goal is to be as fast as the [readfq](https://github.com/lh3/readfq) C library at parsing FASTX files and much (i.e. 25 times) faster than equivalent Python implementations at _k_-mer counting.

For example, a simple Needletail script can count all the bases in a [2.1 gigabyte HiSeq 2500 FASTQ file](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1749083) in 1.1 seconds while a comparable parser with `readfq` takes 2.6 seconds and Biopython takes over one minute (see `bench` folder; measured with `%timeit -r 3 -n 3`, `%timeit -r 3 -n 1` for Biopython). These speed improvements hold for [large FASTQ files](http://www.ebi.ac.uk/ena/data/view/ERX150470) as well.

|                            | needletail  | readfq  | Biopython  |
|----------------------------|---|---|---|
| Mid 2012 MacBook Pro (2GB FASTQ) | 1.83s   | 2.48s  | 2m43s   |
| AWS EC2 r3.xlarge (2GB FASTQ)    | 1.10s  | 2.59s  | 1m47s  |
| AWS EC2 d2.2xlarge (2GB FASTQ)   | 0.93s   | 2.56s  | 1m24s  |
| AWS EC2 d2.2xlarge (55GB FASTQ)   | 34.7s   | 1m6s  | &mdash;  |

_Note: `gcc` with the `-O3` flag was used for `readfq` (`clang -O3` was slower on all tested machines and not used). `rustc` 1.15.1 was used on all machines._

# Example

```rust
extern crate needletail;
use std::env;
use needletail::{fastx};

fn main() {
  let filename: String = env::args().nth(1).unwrap();

  let mut n_bases = 0;
  let mut n_valid_kmers = 0;
  fastx::fastx_file(&filename[..], |seq| {
    // seq.id is the name of the record
    // seq.seq is the base sequence
    // seq.qual is an optional quality score

    // keep track of the total number of bases
    n_bases += seq.seq.len();
    
    // keep track of the number of AAAA (or TTTT via canonicalization) in the 
    /// file (normalize makes sure ever base is capitalized for comparison)
    for (_, kmer, _) in seq.normalize(false).kmers(4, true) {
      if kmer == b"AAAA" {
        n_valid_kmers += 1;
      }
    }
  });
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
needletail = "^0.1.0"
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
