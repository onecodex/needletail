#![crate_name = "needletail"]
//! Needletail is a crate to quickly and easily parse FASTA and FASTQ sequences out of
//! streams/files and manipulate and analyse that data.
//!
//! A contrived example of how to use it:
//! ```
//! extern crate needletail;
//! use needletail::{parse_fastx_file, Sequence, FastxReader};
//!
//! fn main() {
//!     let filename = "tests/data/28S.fasta";
//!
//!     let mut n_bases = 0;
//!     let mut n_valid_kmers = 0;
//!     let mut reader = parse_fastx_file(&filename).expect("valid path/file");
//!     while let Some(record) = reader.next() {
//!         let seqrec = record.expect("invalid record");
//!         // keep track of the total number of bases
//!         n_bases += seqrec.num_bases();
//!         // normalize to make sure all the bases are consistently capitalized and
//!         // that we remove the newlines since this is FASTA
//!         let norm_seq = seqrec.normalize(false);
//!         // we make a reverse complemented copy of the sequence first for
//!         // `canonical_kmers` to draw the complemented sequences from.
//!         let rc = norm_seq.reverse_complement();
//!         // now we keep track of the number of AAAAs (or TTTTs via
//!         // canonicalization) in the file; note we also get the position (i.0;
//!         // in the event there were `N`-containing kmers that were skipped)
//!         // and whether the sequence was complemented (i.2) in addition to
//!         // the canonical kmer (i.1)
//!         for (_, kmer, _) in norm_seq.canonical_kmers(4, &rc) {
//!             if kmer == b"AAAA" {
//!                 n_valid_kmers += 1;
//!             }
//!         }
//!     }
//!     println!("There are {} bases in your file.", n_bases);
//!     println!("There are {} AAAAs in your file.", n_valid_kmers);
//! }
//! ```
//!
#[cfg(any(feature = "python", feature = "python_test"))]
extern crate pyo3;

pub mod bitkmer;
pub mod kmer;
pub mod parser;
pub mod sequence;

pub mod errors;

#[cfg(any(feature = "python", feature = "python_test"))]
pub mod python;

pub use parser::{parse_fastx_file, FastxReader};
pub use sequence::Sequence;
