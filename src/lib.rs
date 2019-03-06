#![crate_name = "needletail"]
extern crate memchr;

#[cfg(feature = "compression")]
extern crate bzip2;
#[cfg(feature = "compression")]
extern crate flate2;
#[cfg(feature = "compression")]
extern crate xz2;
#[cfg(feature = "compression")]
extern crate zip;

pub mod bitkmer;
mod buffer;
pub mod fastx;
pub mod kmer;
pub mod seq;
mod util;
