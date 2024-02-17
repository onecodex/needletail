use std::io;

use memchr::memchr;

use crate::errors::ParseError;
use crate::parser::record::SequenceRecord;

pub(crate) const BUFSIZE: usize = 64 * 1024;

/// Remove a final '\r' from a byte slice
#[inline]
pub(crate) fn trim_cr(line: &[u8]) -> &[u8] {
    if let Some((&b'\r', remaining)) = line.split_last() {
        remaining
    } else {
        line
    }
}

/// Standard buffer policy: buffer size
/// doubles until it reaches 8 MiB. Above, it will
/// increase in steps of 8 MiB. Buffer size is not limited,
/// it could theoretically grow indefinitely.
pub(crate) fn grow_to(current_size: usize) -> usize {
    if current_size < 1 << 23 {
        current_size * 2
    } else {
        current_size + (1 << 23)
    }
}

/// Makes sure the buffer is full after this call (unless EOF reached)
/// code adapted from `io::Read::read_exact`
pub(crate) fn fill_buf<R>(reader: &mut buffer_redux::BufReader<R>) -> io::Result<usize>
where
    R: io::Read,
{
    let initial_size = reader.buffer().len();
    let mut num_read = 0;
    while initial_size + num_read < reader.capacity() {
        match reader.read_into_buf() {
            Ok(0) => break,
            Ok(n) => num_read += n,
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(num_read)
}

/// Holds line number and byte offset of our current state in a parser
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Position {
    pub(crate) line: u64,
    pub(crate) byte: u64,
}

impl Position {
    pub fn new(line: u64, byte: u64) -> Self {
        Self { line, byte }
    }

    /// Line number (starting with 1)
    pub fn line(&self) -> u64 {
        self.line
    }

    /// Byte offset within the file
    pub fn byte(&self) -> u64 {
        self.byte
    }
}

/// FASTA or FASTQ?
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum Format {
    Fasta,
    Fastq,
}

impl Format {
    pub fn start_char(&self) -> char {
        match self {
            Self::Fasta => '>',
            Self::Fastq => '@',
        }
    }
}

/// Whether it uses \r\n or only \n
#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum LineEnding {
    Windows,
    Unix,
}

impl LineEnding {
    pub fn to_bytes(&self) -> Vec<u8> {
        match self {
            Self::Windows => vec![b'\r', b'\n'],
            Self::Unix => vec![b'\n'],
        }
    }
}

pub fn find_line_ending(bytes: &[u8]) -> Option<LineEnding> {
    if !bytes.is_empty() {
        if let Some(idx) = memchr(b'\n', bytes) {
            if idx > 0 && bytes[idx - 1] == b'\r' {
                return Some(LineEnding::Windows);
            }

            return Some(LineEnding::Unix);
        }
    }
    None
}
/// The main trait, iterator-like, that the FASTA and FASTQ readers implement
pub trait FastxReader: Send {
    /// Gets the next record in the stream.
    /// This imitates the Iterator API but does not support any iterator functions.
    /// This returns None once we reached the EOF.
    fn next(&mut self) -> Option<Result<SequenceRecord, ParseError>>;
    /// Returns the current line/byte in the stream we are reading from
    fn position(&self) -> &Position;
    /// Returns whether the current stream uses Windows or Unix style line endings
    /// It is `None` only before calling `next`, once `next` has been called it will always
    /// return a line ending.
    fn line_ending(&self) -> Option<LineEnding>;
}
