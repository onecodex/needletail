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

// mod fasta;
pub mod fastq;

use std::cmp::min;
use std::io::{Cursor, Read};
use std::str;

#[cfg(feature = "compression")]
use bzip2::read::BzDecoder;
#[cfg(feature = "compression")]
use flate2::read::MultiGzDecoder;
#[cfg(feature = "compression")]
use xz2::read::XzDecoder;

use crate::buffer::RecBuffer;
// pub use crate::formats::fasta::FASTA;
pub use crate::formats::fastq::{get_fastq, FASTQ};
use crate::seq::Sequence;
use crate::util::{ParseError, ParseErrorType};

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
    let mut first = vec![0];
    reader.read_exact(&mut first)?;
    let ft = match first[0] {
        b'>' => "FASTA",
        b'@' => "FASTQ",
        _ => "",
    };
    type_callback(ft);

    let mut rec_reader = RecBuffer::new(reader, 10_000_000, &first)?;
    match ft {
        "FASTA" => {
            // while let Some(s) = rec_reader.next::<FASTA>() {
            //     callback(Sequence::from(s?));
            // }
        },
        "FASTQ" => {
            while let Some(s) = get_fastq(&mut rec_reader) {
                callback(Sequence::from(s?));
            }
        },
        _ => {
            return Err(
                ParseError::new("Bad starting byte", ParseErrorType::InvalidHeader)
                    .record(0)
                    .context(String::from_utf8_lossy(&first)),
            )
        },
    };

    if !rec_reader.last {
        return Err(
            ParseError::new("File ended abruptly", ParseErrorType::PrematureEOF)
                .record(rec_reader.count),
        );
    }
    for c in &rec_reader.buf[rec_reader.pos..] {
        if c != &b'\r' && c != &b'\n' {
            let end = min(rec_reader.pos + 16, rec_reader.buf.len());
            let context = String::from_utf8_lossy(&rec_reader.buf[rec_reader.pos..end]);
            return Err(ParseError::new(
                "File had extra data past end of records",
                ParseErrorType::PrematureEOF,
            )
            .record(rec_reader.count)
            .context(context));
        }
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
