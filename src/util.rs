use std::borrow::Cow;
use std::error;
use std::fmt;
use std::io;
use std::str;

use memchr::{memchr, memchr2};

#[cfg(feature = "compression")]
use zip::result::ZipError;

#[derive(Clone, Debug, PartialEq)]
pub enum ParseError {
    PrematureEOF,
    Invalid(String),
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let msg = match *self {
            ParseError::PrematureEOF => "File ended prematurely",
            ParseError::Invalid(ref s) => &s,
        };
        write!(f, "{}", msg)
    }
}

impl error::Error for ParseError {
    fn description(&self) -> &str {
        "ParseError"
    }

    fn cause(&self) -> Option<&error::Error> {
        None
    }
}

impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> ParseError {
        ParseError::Invalid(err.to_string())
    }
}

impl From<str::Utf8Error> for ParseError {
    fn from(err: str::Utf8Error) -> ParseError {
        ParseError::Invalid(err.to_string())
    }
}

#[cfg(feature = "compression")]
impl From<ZipError> for ParseError {
    fn from(err: ZipError) -> ParseError {
        ParseError::Invalid(err.to_string())
    }
}

/// remove newlines from within FASTX records; currently the rate limiting step
/// in FASTX parsing (in general; readfq also exhibits this behavior)
#[inline]
pub fn strip_whitespace<'a>(seq: &'a [u8]) -> Cow<'a, [u8]> {
    let mut new_buf = Vec::with_capacity(seq.len());
    let mut i = 0;
    while i < seq.len() {
        match memchr2(b'\r', b'\n', &seq[i..]) {
            None => {
                new_buf.extend_from_slice(&seq[i..]);
                break;
            },
            Some(match_pos) => {
                new_buf.extend_from_slice(&seq[i..i + match_pos]);
                i += match_pos + 1;
            },
        }
    }
    Cow::Owned(new_buf)
}

/// Like memchr, but handles a two-byte sequence (unlike memchr::memchr2, this
/// looks for the bytes in sequence not either/or).
///
/// Also returns if any other `b1`s were found in the sequence
#[inline]
pub fn memchr_both(b1: u8, b2: u8, seq: &[u8]) -> Option<usize> {
    // TODO: b2 is going to be much rarer for us in FASTAs, so we should search for that instead
    // (but b1 is going to be the rarer character in FASTQs so this is optimized for that and we
    // should allow a choice between the two)
    let mut pos = 0;
    loop {
        match memchr(b1, &seq[pos..]) {
            None => return None,
            Some(match_pos) => {
                if pos + match_pos + 1 == seq.len() {
                    return None;
                } else if seq[pos + match_pos + 1] == b2 {
                    return Some(pos + match_pos);
                } else {
                    pos += match_pos + 1;
                }
            },
        }
    }
}

#[test]
fn test_memchr_both() {
    let pos = memchr_both(b'\n', b'-', &b"test\n-this"[..]);
    assert_eq!(pos, Some(4));

    let pos = memchr_both(b'\n', b'-', &b"te\nst\n-this"[..]);
    assert_eq!(pos, Some(5));
}
