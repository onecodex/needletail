#![crate_name = "needletail"]
pub mod bitkmer;
mod buffer;
pub mod fastx;
pub mod kmer;
pub mod seq;
mod util;

pub use fastx::parse_sequences;
