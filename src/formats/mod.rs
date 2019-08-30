//! This module contains functions for reading FASTA data from either within
//! memory or from files.
//!
//! # Design Note
//!
//! These functions are designed to take callbacks to process the FASTX records
//! they read. It would be nice to present a FASTX Iterator that downstream users
//! can use at some point, but this is only possible for in-memory data (using
//! MemProducer from nom) because otherwise the lifetime of each record is linked
//! to what part of the file we're reading and Rust doesn't support "streaming
//! iterators" otherwise. Maybe some day.
//!
//! See: https://github.com/emk/rust-streaming

mod buffer;
mod fasta;
mod fastq;

use std::cmp::min;
use std::io::{Cursor, Read};
use std::str;

#[cfg(feature = "compression")]
use bzip2::read::BzDecoder;
#[cfg(feature = "compression")]
use flate2::read::MultiGzDecoder;
#[cfg(feature = "compression")]
use xz2::read::XzDecoder;

use crate::formats::buffer::RecReader;
pub use crate::formats::fasta::FASTA;
pub use crate::formats::fastq::FASTQ;
use crate::seq::Sequence;
use crate::util::{ParseError, ParseErrorType};

fn check_end<T>(rec_reader: &RecReader<T>) -> Result<(), ParseError> {
    // check if there's anything left stuff in the buffer (besides returns)
    let rec_buffer = rec_reader.get_buffer();
    if !rec_buffer.last {
        return Err(
            ParseError::new("File ended abruptly", ParseErrorType::PrematureEOF)
                .record(rec_buffer.count),
        );
    }
    for c in &rec_buffer.buf[rec_buffer.pos..] {
        if c != &b'\r' && c != &b'\n' {
            let end = min(rec_buffer.pos + 16, rec_buffer.buf.len());
            let context = String::from_utf8_lossy(&rec_buffer.buf[rec_buffer.pos..end]);
            return Err(ParseError::new(
                "File had extra data past end of records",
                ParseErrorType::PrematureEOF,
            )
            .record(rec_buffer.count + 1)
            .context(context));
        }
    }
    Ok(())
}

/// Internal function abstracting over byte and file FASTX parsing
fn seq_reader<F, R, T>(
    reader: &mut R,
    mut callback: F,
    type_callback: &mut T,
) -> Result<(), ParseError>
where
    F: for<'a> FnMut(Sequence<'a>) -> (),
    R: Read,
    T: ?Sized + FnMut(&'static str) -> (),
{
    // infer the type of the sequencing data
    let mut first = vec![0];
    reader.read_exact(&mut first)?;
    let file_type = match first[0] {
        b'>' => Ok("FASTA"),
        b'@' => Ok("FASTQ"),
        _ => Err(
            ParseError::new("Could not detect file type", ParseErrorType::InvalidHeader)
                .record(0)
                .context(String::from_utf8_lossy(&first)),
        ),
    }?;
    type_callback(file_type);

    match file_type {
        "FASTA" => {
            let mut rec_reader = RecReader::<FASTA>::new(reader, 10_000_000, &first)?;
            loop {
                let used = {
                    let mut rec_buffer = rec_reader.get_buffer();
                    for s in rec_buffer.by_ref() {
                        callback(Sequence::from(s?));
                    }
                    rec_buffer.used()
                };
                if rec_reader.refill(used)? {
                    break;
                }
            }
            check_end(&rec_reader)?;
        },
        "FASTQ" => {
            let mut rec_reader = RecReader::<FASTQ>::new(reader, 10_000_000, &first)?;
            loop {
                let used = {
                    let mut rec_buffer = rec_reader.get_buffer();
                    for s in rec_buffer.by_ref() {
                        callback(Sequence::from(s?));
                    }
                    rec_buffer.used()
                };
                if rec_reader.refill(used)? {
                    break;
                }
            }
            check_end(&rec_reader)?;
        },
        _ => panic!("A file type was inferred that could not be parsed"),
    }
    Ok(())
}

#[cfg(not(feature = "compression"))]
pub fn parse_sequences<F, R, T>(
    mut reader: R,
    mut type_callback: T,
    callback: F,
) -> Result<(), ParseError>
where
    F: for<'a> FnMut(Sequence<'a>) -> (),
    R: Read,
    T: FnMut(&'static str) -> (),
{
    //! Opens a `Read` stream and parses the FASTX records out. Also takes a "type_callback"
    //! that gets called as soon as we determine if the records are FASTA or FASTQ.
    fastx_reader(&mut reader, callback, &mut type_callback)
}

#[cfg(feature = "compression")]
pub fn parse_sequences<F, R, T>(
    mut reader: R,
    mut type_callback: T,
    callback: F,
) -> Result<(), ParseError>
where
    F: for<'a> FnMut(Sequence<'a>) -> (),
    R: Read,
    T: FnMut(&'static str) -> (),
{
    //! Opens a `Read` stream and parses the FASTX records out. Also takes a "type_callback"
    //! that gets called as soon as we determine if the records are FASTA or FASTQ.
    //!  If a file starts with a gzip or other header, transparently decompress it.
    let mut first = vec![0];
    reader.read_exact(&mut first)?;
    if first[0] == 0x1F {
        // gz files
        reader.read_exact(&mut first)?;
        if first[0] != 0x8B {
            return Err(ParseError::new(
                "Bad gz header",
                ParseErrorType::BadCompression,
            ));
        }
        let cursor = Cursor::new(vec![0x1F, 0x8B]);
        let mut gz_reader = MultiGzDecoder::new(cursor.chain(reader));

        seq_reader(&mut gz_reader, callback, &mut type_callback)
    } else if first[0] == 0x42 {
        // bz files
        reader.read_exact(&mut first)?;
        if first[0] != 0x5A {
            return Err(ParseError::new(
                "Bad bz header",
                ParseErrorType::BadCompression,
            ));
        }
        let cursor = Cursor::new(vec![0x42, 0x5A]);
        let mut bz_reader = BzDecoder::new(cursor.chain(reader));

        seq_reader(&mut bz_reader, callback, &mut type_callback)
    } else if first[0] == 0xFD {
        // xz files
        reader.read_exact(&mut first)?;
        if first[0] != 0x37 {
            return Err(ParseError::new(
                "Bad xz header",
                ParseErrorType::BadCompression,
            ));
        }
        let cursor = Cursor::new(vec![0xFD, 0x37]);
        let mut xz_reader = XzDecoder::new(cursor.chain(reader));

        seq_reader(&mut xz_reader, callback, &mut type_callback)
    } else {
        let cursor = Cursor::new(first);
        let mut reader = cursor.chain(reader);
        seq_reader(&mut reader, callback, &mut type_callback)
    }
}
