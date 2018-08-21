#![crate_name = "needletail"]
extern crate memchr;

#[cfg(feature = "compression")]
extern crate flate2;
#[cfg(feature = "compression")]
extern crate bzip2;
#[cfg(feature = "compression")]
extern crate xz2;
#[cfg(feature = "compression")]
extern crate zip;

pub mod fastx;
pub mod kmer;
pub mod bitkmer;
pub mod seq;
mod buffer;
mod util;
