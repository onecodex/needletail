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

use std::io::{Cursor, Read};
use std::str;

#[cfg(feature = "compression")]
use bzip2::read::BzDecoder;
#[cfg(feature = "compression")]
use flate2::read::MultiGzDecoder;
#[cfg(feature = "compression")]
use xz2::read::XzDecoder;

pub use crate::formats::buffer::{RecBuffer, RecParser};
pub use crate::formats::fasta::{FastaRecord, FastaParser};
pub use crate::formats::fastq::{FastqRecord, FastqParser};
use crate::seq::Sequence;
use crate::util::{ParseError, ParseErrorType};

#[macro_export]
macro_rules! parse_stream {
    ($reader:expr, $first:expr, $reader_type: ty, $rec: ident, $handler: block) => {{
        use $crate::formats::{RecBuffer, RecParser};
        let mut buffer = RecBuffer::new($reader, 500_000, &$first)?;
        let mut rec_reader = <$reader_type>::from_buffer(&buffer.buf, buffer.last);
        // TODO: do something with the header?
        let mut record_count: usize = 0;
        rec_reader.header().map_err(|e| e.record(record_count))?;
        let used = rec_reader.used();
        if !buffer.refill(used).map_err(|e| e.record(record_count))? {
            loop {
                let used = {
                    let mut rec_reader = <$reader_type>::from_buffer(&buffer.buf, buffer.last);
                    for s in rec_reader.by_ref() {
                        record_count += 1;
                        let $rec = s.map_err(|e| e.record(record_count))?;
                        $handler
                    }
                    rec_reader.used()
                };
                if buffer.refill(used).map_err(|e| e.record(record_count))? {
                    break;
                }
            }
        }
        let rec_reader = <$reader_type>::from_buffer(&buffer.buf, buffer.last);
        rec_reader.eof().map_err(|e| e.record(record_count + 1))?;
    }};
}

/// Internal function abstracting over byte and file FASTX parsing
#[inline]
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
        "FASTA" => parse_stream!(reader, first, FastaParser, rec, {
            callback(Sequence::from(rec))
        }),
        "FASTQ" => parse_stream!(reader, first, FastqParser, rec, {
            callback(Sequence::from(rec))
        }),
        _ => panic!("A file type was inferred that could not be parsed"),
    };
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
    //! If a file starts with a gzip or other header, transparently decompress it.
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
