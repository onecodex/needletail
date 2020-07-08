//! Handles all the FASTA/FASTQ parsing
use std::fs::File;
use std::io::{stdin, Cursor, Read};
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

// Magic bytes for each compression format
#[cfg(feature = "compression")]
const GZ_MAGIC: [u8; 2] = [0x1F, 0x8B];
#[cfg(feature = "compression")]
const BZ_MAGIC: [u8; 2] = [0x42, 0x5A];
#[cfg(feature = "compression")]
const XZ_MAGIC: [u8; 2] = [0xFD, 0x37];

fn get_fastx_reader<'a, R: 'a + io::Read + Send>(
    reader: R,
    first_byte: u8,
) -> Result<Box<dyn FastxReader + 'a>, ParseError> {
    match first_byte {
        b'>' => Ok(Box::new(FastaReader::new(reader))),
        b'@' => Ok(Box::new(FastqReader::new(reader))),
        _ => Err(ParseError::new_unknown_format(first_byte)),
    }
}

/// The main entry point of needletail if you're reading from something that impls std::io::Read
/// This automatically detects whether the file is:
/// 1. compressed: gzip, bz and xz are supported and will use the appropriate decoder
/// 2. FASTA or FASTQ: the right parser will be automatically instantiated
/// 1 is only available if the `compression` feature is enabled.
pub fn parse_fastx_reader<'a, R: 'a + io::Read + Send>(
    mut reader: R,
) -> Result<Box<dyn FastxReader + 'a>, ParseError> {
    let mut first_two_bytes = [0; 2];
    reader.read_exact(&mut first_two_bytes)?;
    let first_two_cursor = Cursor::new(first_two_bytes);
    let new_reader = first_two_cursor.chain(reader);

    match first_two_bytes {
        #[cfg(feature = "compression")]
        GZ_MAGIC => {
            let mut gz_reader = MultiGzDecoder::new(new_reader);
            let mut first = [0; 1];
            gz_reader.read_exact(&mut first)?;
            let r = Cursor::new(first).chain(gz_reader);
            get_fastx_reader(r, first[0])
        }
        #[cfg(feature = "compression")]
        BZ_MAGIC => {
            let mut bz_reader = BzDecoder::new(new_reader);
            let mut first = [0; 1];
            bz_reader.read_exact(&mut first)?;
            let r = Cursor::new(first).chain(bz_reader);
            get_fastx_reader(r, first[0])
        }
        #[cfg(feature = "compression")]
        XZ_MAGIC => {
            let mut xz_reader = XzDecoder::new(new_reader);
            let mut first = [0; 1];
            xz_reader.read_exact(&mut first)?;
            let r = Cursor::new(first).chain(xz_reader);
            get_fastx_reader(r, first[0])
        }
        _ => get_fastx_reader(new_reader, first_two_bytes[0]),
    }
}

/// The main entry point of needletail if you're reading from stdin.
/// Shortcut to calling `parse_fastx_reader` with `stdin()`
pub fn parse_fastx_stdin() -> Result<Box<dyn FastxReader>, ParseError> {
    let stdin = stdin();
    parse_fastx_reader(stdin)
}

/// The main entry point of needletail if you're reading from a file.
/// Shortcut to calling `parse_fastx_reader` with a file
pub fn parse_fastx_file<P: AsRef<Path>>(path: P) -> Result<Box<dyn FastxReader>, ParseError> {
    parse_fastx_reader(File::open(&path)?)
}

pub use record::{mask_header_tabs, mask_header_utf8, write_fasta, write_fastq, SequenceRecord};
use std::io;
pub use utils::{Format, LineEnding};
