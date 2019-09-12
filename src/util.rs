use std::error;
use std::fmt;
use std::io;
use std::str;

use memchr::memchr_iter;

/// The type of error that occured during file parsing
#[derive(Clone, Debug, PartialEq)]
pub enum ParseErrorType {
    /// An error was encountered parsing a record's header
    InvalidHeader,
    /// An error was encountered parsing a record
    InvalidRecord,
    /// An error happened during file/stream input/output
    IOError,
    /// A generic error occured
    Invalid,
}

/// The only error type that needletail returns
#[derive(Clone, Debug, PartialEq)]
pub struct ParseError {
    /// The number of the record where the error occured (if relevant)
    pub record: usize,
    /// Grep-able context for where the error occured (if relevant)
    pub context: String,
    /// A description of what went wrong
    pub msg: String,
    /// The type of error that occured
    pub error_type: ParseErrorType,
}

impl ParseError {
    pub fn new<S>(msg: S, error_type: ParseErrorType) -> Self
    where
        S: Into<String>,
    {
        ParseError {
            record: 0,
            context: "".to_string(),
            msg: msg.into(),
            error_type,
        }
    }

    pub fn record(mut self, record_number: usize) -> Self {
        self.record = record_number;
        self
    }

    pub fn context<S>(mut self, context: S) -> Self
    where
        S: ToString,
    {
        self.context = context.to_string();
        self
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let msg = match self.error_type {
            ParseErrorType::InvalidHeader => "Invalid record header",
            ParseErrorType::InvalidRecord => "Invalid record content",
            ParseErrorType::IOError => "I/O Error",
            ParseErrorType::Invalid => "",
        };
        write!(f, "{}: {}", msg, self.msg)
    }
}

impl error::Error for ParseError {
    fn description(&self) -> &str {
        "ParseError"
    }

    fn cause(&self) -> Option<&dyn error::Error> {
        None
    }
}

impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> ParseError {
        ParseError::new(err.to_string(), ParseErrorType::IOError)
    }
}

impl From<str::Utf8Error> for ParseError {
    fn from(err: str::Utf8Error) -> ParseError {
        ParseError::new(err.to_string(), ParseErrorType::Invalid)
    }
}

/// Like memchr, but handles a two-byte sequence (unlike memchr::memchr2, this
/// looks for the bytes in sequence not either/or).
#[inline]
pub fn memchr_both(b1: u8, b2: u8, seq: &[u8]) -> Option<usize> {
    for idx in memchr_iter(b1, &seq) {
        if idx + 1 < seq.len() && seq[idx + 1] == b2 {
            return Some(idx);
        }
    }
    None
}

/// Like `memchr_both`, but searches for `b2` first and then checks if b1 is also
/// there. Should give the same results as `memchr_both` but will be faster if
/// `b2` is a rarer character than `b1` (such as `>` vs `\n` in multiline FASTAs).
#[inline]
pub fn memchr_both_last(b1: u8, b2: u8, seq: &[u8]) -> Option<usize> {
    for idx in memchr_iter(b2, &seq) {
        if idx != 0 && seq[idx - 1] == b1 {
            return Some(idx - 1);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_memchr_both() {
        let pos = memchr_both(b'\n', b'-', &b"test\n-this"[..]);
        assert_eq!(pos, Some(4));

        let pos = memchr_both(b'\n', b'-', &b"te\nst\n-this"[..]);
        assert_eq!(pos, Some(5));
    }

    #[test]
    fn test_memchr_both_last() {
        let pos = memchr_both_last(b'\n', b'-', &b"test\n-this"[..]);
        assert_eq!(pos, Some(4));

        let pos = memchr_both_last(b'\n', b'-', &b"te\nst\n-this"[..]);
        assert_eq!(pos, Some(5));

        let pos = memchr_both_last(b'\n', b'-', &b"-te\nst\n-this"[..]);
        assert_eq!(pos, Some(6));
    }

}
