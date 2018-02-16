use std::borrow::Cow;
use std::error;
use std::fmt;
use std::io;
use std::str;

use memchr::{memchr, memchr2};


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




/// remove newlines from within FASTX records; currently the rate limiting step
/// in FASTX parsing (in general; readfq also exhibits this behavior)
#[inline]
fn strip_whitespace<'a>(seq: &'a [u8], has_newlines: bool) -> Cow<'a, [u8]> {
    // if we already know there are no newlines inside, we only have to remove
    // ones at the end
    if !has_newlines {
        let mut s_seq = seq;
        if s_seq[s_seq.len() - 1] == b'\n' {
            s_seq = &s_seq[..s_seq.len() - 1];
        }
        if s_seq[s_seq.len() - 1] == b'\r' {
            s_seq = &s_seq[..s_seq.len() - 1];
        }
        return Cow::Borrowed(s_seq);
    }

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
pub fn memchr_both(b1: u8, b2: u8, seq: &[u8]) -> (Option<usize>, bool) {
    let mut pos = 0;
    let mut found_newline = false;
    loop {
        match memchr(b1, &seq[pos..]) {
            None => return (None, found_newline),
            Some(match_pos) => {
                if pos + match_pos + 1 == seq.len() {
                    return (None, found_newline);
                } else if seq[pos + match_pos + 1] == b2 {
                    return (Some(pos + match_pos), found_newline);
                } else {
                    pos += match_pos + 1;
                    found_newline = true;
                }
            },
        }
    }
}


#[test]
fn test_memchr_both() {
    let (pos, newline) = memchr_both(b'\n', b'-', &b"test\n-this"[..]);
    assert_eq!(pos, Some(4));
    assert_eq!(newline, false);

    let (pos, newline) = memchr_both(b'\n', b'-', &b"te\nst\n-this"[..]);
    assert_eq!(pos, Some(5));
    assert_eq!(newline, true);
}
