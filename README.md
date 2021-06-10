![CI](https://github.com/onecodex/needletail/workflows/CI/badge.svg)
[![crates.io](https://img.shields.io/crates/v/needletail.svg)](https://crates.io/crates/needletail)

# Needletail

Needletail is a MIT-licensed, minimal-copying FASTA/FASTQ parser and _k_-mer processing library for Rust.

The goal is to write a fast *and* well-tested set of functions that more specialized bioinformatics programs can use.
Needletail's goal is to be as fast as the [readfq](https://github.com/lh3/readfq) C library at parsing FASTX files and much (i.e. 25 times) faster than equivalent Python implementations at _k_-mer counting.

## Example

```rust
extern crate needletail;
use needletail::{parse_fastx_file, Sequence, FastxReader};

fn main() {
    let filename = "tests/data/28S.fasta";

    let mut n_bases = 0;
    let mut n_valid_kmers = 0;
    let mut reader = parse_fastx_file(&filename).expect("valid path/file");
    while let Some(record) = reader.next() {
        let seqrec = record.expect("invalid record");
        // keep track of the total number of bases
        n_bases += seqrec.num_bases();
        // normalize to make sure all the bases are consistently capitalized and
        // that we remove the newlines since this is FASTA
        let norm_seq = seqrec.normalize(false);
        // we make a reverse complemented copy of the sequence first for
        // `canonical_kmers` to draw the complemented sequences from.
        let rc = norm_seq.reverse_complement();
        // now we keep track of the number of AAAAs (or TTTTs via
        // canonicalization) in the file; note we also get the position (i.0;
        // in the event there were `N`-containing kmers that were skipped)
        // and whether the sequence was complemented (i.2) in addition to
        // the canonical kmer (i.1)
        for (_, kmer, _) in norm_seq.canonical_kmers(4, &rc) {
            if kmer == b"AAAA" {
                n_valid_kmers += 1;
            }
        }
    }
    println!("There are {} bases in your file.", n_bases);
    println!("There are {} AAAAs in your file.", n_valid_kmers);
}
```

## Installation

Needletail requires `rust` and `cargo` to be installed.
Please use either your local package manager (`homebrew`, `apt-get`, `pacman`, etc) or install these via [rustup](https://www.rustup.rs/).

Once you have Rust set up, you can include needletail in your `Cargo.toml` file like:
```shell
[dependencies]
needletail = "0.4"
```

To install needletail itself for development:
```shell
git clone https://github.com/onecodex/needletail
cargo test  # to run tests
```

### Python

#### Documentation

For a real example, you can refer to `test_python.py`.

The python library only raise one type of exception: `NeedletailError`.

There are 2 ways to parse a FASTA/FASTQ: one if you have a string (`parse_fastx_string(content: str)`) or a path to a file
(`parse_fastx_file(path: str)`). Those functions will raise if the file is not found or if the content is invalid and will return
an iterator.


```python
from needletail import parse_fastx_file, NeedletailError, reverse_complement, normalize_seq

try:
    for record in parse_fastx_file("myfile.fastq"):
        print(record.id)
        print(record.seq)
        print(record.qual)
except NeedletailError:
    print("Invalid Fastq file")
```

A record has the following shape:

```python
class Record:
    id: str
    seq: str
    qual: Optional[str]

    def is_fasta(self) -> bool
    def is_fastq(self) -> bool
    def normalize(self, iupac: bool)
```

Note that `normalize` (see <https://docs.rs/needletail/0.4.1/needletail/sequence/fn.normalize.html> for what it does) will mutate `self.seq`.
It is also available as the `normalize_seq(seq: str, iupac: bool)` function which will return the normalized sequence in this case.

Lastly, there is also a `reverse_complement(seq: str)` that will do exactly what it says. This will not raise an error if you pass some invalid
characters.

#### Building

To work on the Python library on a Mac OS X/Unix system (requires Python 3):
```bash
pip install maturin

# finally, install the library in the local virtualenv
maturin develop --cargo-extra-args="--features=python"
```

To build the binary wheels and pushing to PyPI

```
# The Mac build requires switching through a few different python versions
maturin build --cargo-extra-args="--features=python" --release --strip

# The linux build is automated through cross-compiling in a docker image
docker run --rm -v $(pwd):/io konstin2/maturin:master build --cargo-extra-args="--features=python" --release --strip
twine upload target/wheels/*
```

## Getting Help

Questions are best directed as GitHub issues. We plan to add more documentation soon, but in the meantime "doc" comments are included in the source.

## Contributing

Please do! We're happy to discuss possible additions and/or accept pull requests.

## Acknowledgements
Starting from 0.4, the parsers algorithms is taken from [seq_io](https://github.com/markschl/seq_io). While it has been slightly modified, it is mainly
coming from that library. Links to the original files are available in `src/parser/fast{a,q}.rs`.
