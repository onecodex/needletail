//! The errors needletail can return; only when parsing FASTA/FASTQ files

use crate::parser::Format;
use std::error::Error as StdError;
use std::{fmt, num};
use std::io;
use std::path::Display;

/// Represents where we were in a file when an error occurred.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ErrorPosition {
    /// Line number where the error occurred (starting with 1)
    pub line: u64,
    /// ID of record if available
    pub id: Option<String>,
}

impl fmt::Display for ErrorPosition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(id) = self.id.as_ref() {
            write!(f, "record '{}' at ", id)?;
        }
        write!(f, "line {}", self.line)
    }
}

/// The type of error that occurred during file parsing
#[derive(Clone, Debug, PartialEq)]
pub enum ParseErrorKind {
    /// An error happened during file/stream input/output
    Io,
    /// The file didn't start with `@` or `>` and we didn't know what to expect yet
    UnknownFormat,
    /// Invalid start byte of record encountered (expected `@` in FASTQ and `>` in FASTA)
    InvalidStart,
    /// The separator line in a FASTQ file is not valid (no `+`)
    InvalidSeparator,
    /// Sequence and quality lengths are not equal (in a FASTQ file only)
    UnequalLengths,
    /// Truncated record found
    UnexpectedEnd,
    /// The file appears to be empty
    EmptyFile,
}

/// The error type of Index/IndexReader operations
#[derive(Clone, Debug, PartialEq)]
pub enum IndexErrorKind {
    /// The fai file IO error
    FaiIo,
    /// The fai file is not correct format
    FaiFormatError,
    /// Unknown sequence name in fai file
    UnknownSeqName,
    /// Region is not valid
    InvalidRegion,
    /// Seek IO error
    Io,
    /// Inner Reader error
    InnerReader,
}

/// The error type of Index/IndexReader operations
#[derive(Clone, Debug, PartialEq)]
pub struct IndexError {
    /// A description of what went wrong
    pub msg: String,
    /// The type of error that occurred
    pub kind: IndexErrorKind,
}

impl IndexError {

    pub fn new_fai_format_err() -> Self {
        IndexError {
            msg: String::from("Fai index format parse error, please check the format."),
            kind:IndexErrorKind::FaiFormatError,
        }
    }

    pub fn new_fai_io_err(path: Display) -> Self {
        IndexError {
            msg: format!("Fai index file of `{}` read error.", path),
            kind:IndexErrorKind::FaiIo,
        }
    }

    pub fn new_seq_name_err(seq_name: &str) -> Self {
        IndexError {
            msg: format!("Unknown sequence name `{}` in fasta file.", seq_name),
            kind:IndexErrorKind::UnknownSeqName,
        }
    }

    pub fn new_invalid_region_err() -> Self {
        IndexError {
            msg: String::from("Invalid query region."),
            kind:IndexErrorKind::InvalidRegion,
        }
    }

    pub fn new_io_err() -> Self {
        IndexError {
            msg: String::from("IO error."),
            kind:IndexErrorKind::Io,
        }
    }
}

/// The only error type that needletail returns
#[derive(Clone, Debug, PartialEq)]
pub struct ParseError {
    /// A description of what went wrong
    pub msg: String,
    /// The type of error that occurred
    pub kind: ParseErrorKind,
    /// Position within file
    pub position: ErrorPosition,
    /// The format of the file we were parsing
    pub format: Option<Format>,
}

impl ParseError {
    pub fn new_invalid_start(byte_found: u8, position: ErrorPosition, format: Format) -> Self {
        let msg = format!(
            "Expected '{}' but found '{}",
            format.start_char(),
            (byte_found as char).escape_default()
        );
        ParseError {
            kind: ParseErrorKind::InvalidStart,
            msg,
            position,
            format: Some(format),
        }
    }

    pub fn new_invalid_separator(byte_found: u8, position: ErrorPosition) -> Self {
        let msg = format!(
            "Expected '+' separator but found '{}",
            (byte_found as char).escape_default()
        );
        ParseError {
            kind: ParseErrorKind::InvalidSeparator,
            msg,
            position,
            format: Some(Format::Fastq),
        }
    }

    pub fn new_unknown_format(byte_found: u8) -> Self {
        let msg = format!(
            "Expected '@' or '>' at the start of the file but found '{}'.",
            (byte_found as char).escape_default()
        );
        ParseError {
            kind: ParseErrorKind::UnknownFormat,
            msg,
            position: ErrorPosition::default(),
            format: Some(Format::Fastq),
        }
    }

    pub fn new_unequal_length(seq_len: usize, qual_len: usize, position: ErrorPosition) -> Self {
        let msg = format!(
            "Sequence length is {} but quality length is {}",
            seq_len, qual_len
        );
        ParseError {
            kind: ParseErrorKind::UnequalLengths,
            msg,
            position,
            format: Some(Format::Fastq),
        }
    }

    pub fn new_unexpected_end(position: ErrorPosition, format: Format) -> Self {
        ParseError {
            msg: String::new(),
            kind: ParseErrorKind::UnexpectedEnd,
            position,
            format: Some(format),
        }
    }

    pub fn new_empty_file() -> Self {
        ParseError {
            msg: String::from("Failed to read the first two bytes. Is the file empty?"),
            kind: ParseErrorKind::EmptyFile,
            position: ErrorPosition::default(),
            format: None,
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.kind {
            ParseErrorKind::Io => write!(f, "I/O error: {}", self.msg),
            ParseErrorKind::UnequalLengths
            | ParseErrorKind::InvalidStart
            | ParseErrorKind::UnknownFormat
            | ParseErrorKind::EmptyFile
            | ParseErrorKind::InvalidSeparator => write!(f, "{} ({})", self.msg, self.position),
            ParseErrorKind::UnexpectedEnd => {
                write!(f, "Unexpected end of input ({}).", self.position)
            }
        }
    }
}

impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> ParseError {
        ParseError {
            msg: err.to_string(),
            kind: ParseErrorKind::Io,
            position: ErrorPosition::default(),
            format: None,
        }
    }
}


impl From<num::ParseIntError> for IndexError {
    fn from(err: num::ParseIntError) -> IndexError {
        IndexError {
            msg: err.to_string(),
            kind: IndexErrorKind::FaiFormatError,
        }
    }
}

impl From<io::Error> for IndexError {
    fn from(err: io::Error) -> IndexError {
        IndexError {
            msg: err.to_string(),
            kind: IndexErrorKind::FaiIo,
        }
    }
}

impl From<ParseError> for IndexError {
    fn from(err: ParseError) -> Self {
        IndexError {
            msg: err.msg,
            kind: IndexErrorKind::InnerReader,
        }
    }
}

impl StdError for ParseError {
    fn cause(&self) -> Option<&dyn StdError> {
        // Ideally we would pass the io::Error but we don't for simplicity sake
        // since we wouldn't be able to `==` on the error kind otherwise.
        // TODO: impl partialeq manually?
        None
    }
}
