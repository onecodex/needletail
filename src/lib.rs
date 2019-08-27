#![crate_name = "needletail"]
pub mod bitkmer;
mod buffer;
pub mod formats;
pub mod kmer;
pub mod seq;
mod util;

pub use formats::parse_sequences;
