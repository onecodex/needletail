//! Handles all the FASTA/FASTQ parsing
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

#[cfg(feature = "compression")]
use bzip2::read::BzDecoder;
#[cfg(feature = "compression")]
use flate2::read::MultiGzDecoder;
#[cfg(feature = "compression")]
use xz2::read::XzDecoder;

use crate::errors::ParseError;
pub use crate::parser::fasta::Reader as FastaReader;
pub use crate::parser::fastq::Reader as FastqReader;

mod record;
mod utils;

mod fasta;
mod fastq;

pub use crate::parser::utils::FastxReader;

// TODO: for now this drops support for stdin parsing. It could be added back if needed though
// by adding a method to the readers to set some initial buffer.

fn read_file(mut f: File) -> Result<Box<dyn FastxReader>, ParseError> {
    let mut first = [0; 1];
    f.read_exact(&mut first)?;
    // Back to the beginning of the file
    f.seek(SeekFrom::Start(0))?;

    match first[0] {
        b'>' => Ok(Box::new(FastaReader::new(f))),
        b'@' => Ok(Box::new(FastqReader::new(f))),
        _ => Err(ParseError::new_unknown_format(first[0])),
    }
}

/// The main entry point of needletail.
/// Parses the file given a path and return an iterator-like reader struct.
/// This automatically detects whether the file is:
/// 1. compressed: gzip, bz and xz are supported and will use the appropriate decoder
/// 2. FASTA or FASTQ: the right parser will be automatically instantiated
/// 1 is only available if the `compression` feature is enabled.
#[cfg(not(feature = "compression"))]
pub fn parse_fastx_file<P: AsRef<Path>>(path: P) -> Result<Box<dyn FastxReader>, ParseError> {
    let f = File::open(&path)?;
    read_file(f)
}

// Magic bytes for each compression format
#[cfg(feature = "compression")]
const GZ_MAGIC: [u8; 2] = [0x1F, 0x8B];
#[cfg(feature = "compression")]
const BZ_MAGIC: [u8; 2] = [0x42, 0x5A];
#[cfg(feature = "compression")]
const XZ_MAGIC: [u8; 2] = [0xFD, 0x37];

/// The main entry point of needletail.
/// Parses the file given a path and return an iterator-like reader struct.
/// This automatically detects whether the file is:
/// 1. compressed: gzip, bz and xz are supported and will use the appropriate decoder
/// 2. FASTA or FASTQ: the right parser will be automatically instantiated
/// 1 is only available if the `compression` feature is enabled.
#[cfg(feature = "compression")]
pub fn parse_fastx_file<P: AsRef<Path>>(path: P) -> Result<Box<dyn FastxReader>, ParseError> {
    let mut f = File::open(&path)?;
    let mut first = [0; 2];
    f.read_exact(&mut first)?;
    // Back to the beginning of the file
    f.seek(SeekFrom::Start(0))?;

    // Not great to say the least, we load the decoder twice since they don't impl Seek.
    // That would cause an issue for compressed gzip and is probably not great for performance
    // Not much compared to decompression overall but still
    match first {
        GZ_MAGIC => {
            let mut gz_reader = MultiGzDecoder::new(f);
            let mut first = [0; 1];
            gz_reader.read_exact(&mut first)?;
            match first[0] {
                b'>' => Ok(Box::new(FastaReader::new(MultiGzDecoder::new(File::open(
                    &path,
                )?)))),
                b'@' => Ok(Box::new(FastqReader::new(MultiGzDecoder::new(File::open(
                    &path,
                )?)))),
                _ => Err(ParseError::new_unknown_format(first[0])),
            }
        }
        BZ_MAGIC => {
            let mut bz_reader = BzDecoder::new(f);
            let mut first = [0; 1];
            bz_reader.read_exact(&mut first)?;
            match first[0] {
                b'>' => Ok(Box::new(FastaReader::new(BzDecoder::new(File::open(
                    &path,
                )?)))),
                b'@' => Ok(Box::new(FastqReader::new(BzDecoder::new(File::open(
                    &path,
                )?)))),
                _ => Err(ParseError::new_unknown_format(first[0])),
            }
        }
        XZ_MAGIC => {
            let mut xz_reader = XzDecoder::new(f);
            let mut first = [0; 1];
            xz_reader.read_exact(&mut first)?;
            match first[0] {
                b'>' => Ok(Box::new(FastaReader::new(XzDecoder::new(File::open(
                    &path,
                )?)))),
                b'@' => Ok(Box::new(FastqReader::new(XzDecoder::new(File::open(
                    &path,
                )?)))),
                _ => Err(ParseError::new_unknown_format(first[0])),
            }
        }
        _ => read_file(f),
    }
}

// TODO: add a method that parses a string but handles decompression as well?

pub use record::{mask_header_tabs, mask_header_utf8, write_fasta, write_fastq, SequenceRecord};
pub use utils::{Format, LineEnding};
