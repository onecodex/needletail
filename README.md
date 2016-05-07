[![Circle CI](https://circleci.com/gh/onecodex/needletail.svg?style=shield&circle-token=65c2b7d87452dba5e8e3e967133311af478632a4)](https://circleci.com/gh/onecodex/needletail)

# Needletail

Needletail is a MIT-licensed, minimal-copying FASTA/FASTQ parser and k-mer processing library.

The goal is to write a fast *and* well-tested set of functions that more-specialized bioinformatics programs can use.
Needletail's goal is to be as fast as the `readfq` C library at parsing FASTX files and much (i.e. 25 times) faster than equivalent Python implementations at k-mer counting.

For example, a simple Needletail script can count all the bases in a 2.1 gigabyte FASTQ file in 2.75 seconds while a comparable parser with `readfq` takes 3.52 seconds (see `bench` folder; measured with `%timeit -r 3 -n 3` on an Early 2015 MacBook Pro).

# Example

```rust
extern crate needletail;
ues std::env;
use needletail::{fastx};

fn main() {
  let filename: String = env::args().nth(1).unwrap();

  let mut n_bases = 0;
  fastx::fastx_file(&filename[..], |seq| {
    // seq.0 is the name of the record
    // seq.1 is the base sequence
    n_bases += seq.1.len();
    // seq.2 is an optional quality score
  });
  println!("There are {} bases in your file.", n_bases);
}
```

# Installation

Needletail requires `rust` and `cargo` to be installed.
Please use either your local package manager (`homebrew`, `apt-get`, `pacman`, etc) or install these via [rustup](https://www.rustup.rs/).

To install needletail itself for development:
```shell
git clone https://github.com/bovee/needletail
cargo test  # to run tests
```

# Getting Help

Questions are best directed as GitHub issues.

Hopefully I'll compile the documentation and put it up as a webpage soon too.

# Contributing

Please do! We're happy to discuss/mentor possible additions and/or accept pull requests.
