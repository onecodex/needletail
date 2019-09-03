#![crate_name = "needletail"]
pub mod bitkmer;
pub mod formats;
pub mod kmer;
pub mod seq;
mod util;

pub use formats::parse_sequences;
pub use util::ParseError;
