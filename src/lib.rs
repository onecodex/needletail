#![crate_name = "needletail"]
extern crate memchr;

#[cfg(feature = "gz")]
extern crate flate2;

pub mod fastx;
pub mod kmer;
pub mod bitkmer;
pub mod seq;
mod buffer;
mod util;
