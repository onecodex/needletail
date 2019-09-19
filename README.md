[![Circle CI](https://circleci.com/gh/onecodex/needletail.svg?style=shield&circle-token=65c2b7d87452dba5e8e3e967133311af478632a4)](https://circleci.com/gh/onecodex/needletail)
[![crates.io](https://img.shields.io/crates/v/needletail.svg)](https://crates.io/crates/needletail)

# Needletail

Needletail is a MIT-licensed, minimal-copying FASTA/FASTQ parser and _k_-mer processing library for Rust.

The goal is to write a fast *and* well-tested set of functions that more specialized bioinformatics programs can use.
Needletail's goal is to be as fast as the [readfq](https://github.com/lh3/readfq) C library at parsing FASTX files and much (i.e. 25 times) faster than equivalent Python implementations at _k_-mer counting.

# Example

```rust
extern crate needletail;
use needletail::{parse_sequence_path, Sequence};
use std::env;

fn main() {
    let filename: String = env::args().nth(1).unwrap();

    let mut n_bases = 0;
    let mut n_valid_kmers = 0;
    parse_sequence_path(
	filename,
        |_| {},
        |seq| {
            // seq.id is the name of the record
            // seq.seq is the base sequence
            // seq.qual is an optional quality score

            // keep track of the total number of bases
            n_bases += seq.seq.len();

            // normalize to make sure all the bases are consistantly capitalized
            let norm_seq = seq.normalize(false);
            // we make a reverse complemented copy of the sequence first for
            // `canonical_kmers` to draw the complemented sequences from.
            let rc = norm_seq.reverse_complement();
            // now we keep track of the number of AAAAs (or TTTTs via
            // canonicalization) in the file; note we also get the postion (i.0;
            // in the event there were `N`-containing kmers that were skipped)
            // and whether the sequence was complemented (i.2) in addition to
            // the canonical kmer (i.1)
            for (_, kmer, _) in norm_seq.canonical_kmers(4, &rc) {
                if kmer == b"AAAA" {
                    n_valid_kmers += 1;
                }
            }
        },
    )
    .expect("parsing failed");
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
needletail = "^0.3.1"
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
