//! Handles all the FASTA/FASTQ parsing
use std::fs::File;
use std::io::{stdin, Cursor, Read};
use std::path::Path;

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

/// The main entry point of needletail.
/// Parses the file given a path and return an iterator-like reader struct.
/// This automatically detects whether the file is:
/// 1. compressed: gzip, bz and xz are supported and will use the appropriate decoder
/// 2. FASTA or FASTQ: the right parser will be automatically instantiated
/// 1 is only available if the `compression` feature is enabled.
pub fn parse_fastx_file<P: AsRef<Path>>(path: P) -> Result<Box<dyn FastxReader>, ParseError> {
    let filename = path.as_ref();
    let f: Box<dyn Read> = if filename == Path::new("-") {
        Box::new(stdin())
    } else {
        Box::new(File::open(path)?)
    };
    parse_fastx_reader(f)
}

/// Parses any input supporting the Read trait and return an iterator-like reader struct.
/// This automatically detects whether the input is:
/// 1. compressed: gzip, bz and xz are supported and will use the appropriate decoder
/// 2. FASTA or FASTQ: the right parser will be automatically instantiated
/// 1 is only available if the `compression` feature is enabled.
pub fn parse_fastx_reader<'a, R: Read + 'a>(
    rdr: R,
) -> Result<Box<dyn FastxReader + 'a>, ParseError> {
    let (mut reader, _) = niffler::get_reader(Box::new(rdr))?;

    let mut first = [0; 1];
    reader.read_exact(&mut first)?;

    let cursor = Cursor::new(first);
    match first[0] {
        b'>' => Ok(Box::new(FastaReader::new(cursor.chain(reader)))),
        b'@' => Ok(Box::new(FastqReader::new(cursor.chain(reader)))),
        _ => Err(ParseError::new_unknown_format(first[0])),
    }
}
// TODO: add a method that parses a string but handles decompression as well?

pub use record::{mask_header_tabs, mask_header_utf8, write_fasta, write_fastq, SequenceRecord};
pub use utils::{Format, LineEnding};
