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

use std::borrow::Cow;
use std::cmp::min;
use std::fs::File;
use std::io::{Cursor, Read, Seek, SeekFrom};
use std::mem;
use std::path::Path;
use std::str;

use memchr::{memchr, memchr2};

use buffer::{RecBuffer, ParseError};

#[cfg(feature = "gz")]
use flate2::read::GzDecoder;

/// A generic FASTX record containing:
///   seq.0 - an id
///   seq.1 - the sequence itself
///   seq.2 - an optional set of quality scores for the sequence
pub type SeqRecord<'a> = (&'a str, Cow<'a, [u8]>, Option<&'a [u8]>);

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
/// Also returns if any other newlines were found in the sequence
#[inline]
fn memchr_both(b1: u8, b2: u8, seq: &[u8]) -> (Option<usize>, bool) {
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


/// Given a RecBuffer, pull out one FASTA record
fn fasta_record<'a>(rb: &'a mut RecBuffer) -> Result<SeqRecord<'a>, ParseError> {
    let _ = rb.mark_field(|buf: &[u8], eof: bool| {
        if buf[0] != b'>' {
            return Err(ParseError::Invalid(String::from("Bad FASTA record")));
        }
        match memchr(b'\n', &buf) {
            None => Err(ParseError::NeedMore),
            Some(pos_end) => Ok(pos_end + 1),
        }
    })?;

    let mut newlines = false;
    let _ = rb.mark_field(|buf: &[u8], eof: bool| {
        match memchr_both(b'\n', b'>', &buf) {
            (None, has_newline) => match eof {
                false => Err(ParseError::NeedMore),
                true => {
                    newlines = has_newline;
                    Ok(buf.len())
                },
            },
            (Some(pos_end), has_newline) => {
                newlines = has_newline;
                Ok(pos_end + 1)
            }
        }
    })?;

    let mut fields: [&[u8]; 2];
    unsafe {
        fields = mem::uninitialized();
    }
    rb.fields(&mut fields);

    let id = match str::from_utf8(&fields[0][1..&fields[0].len() - 1]) { 
        Ok(v) => Ok(v),
        Err(_) => Err(ParseError::Invalid(String::from("FASTA header not UTF8"))),
    }?;

    Ok((id, strip_whitespace(fields[1], newlines), None))
}


/// Given a RecBuffer, pull out one FASTQ record
fn fastq_record<'a>(rb: &'a mut RecBuffer) -> Result<SeqRecord<'a>, ParseError> {
    let _ = rb.mark_field(|buf: &[u8], eof: bool| {
        if buf[0] != b'@' {
            // we allow up to 4 whitespace characters at the very end
            // b/c some FASTQs have extra returns
            if buf.len() <= 4 {
                if !eof {
                    // try to get an eof
                    return Err(ParseError::NeedMore);
                }
                if buf.iter().all(|&c| c == b'\r' || c == b'\n') {
                    return Err(ParseError::EOF);
                }
            }
            return Err(ParseError::Invalid(String::from("Bad FASTQ record")));
        }
        match memchr(b'\n', &buf) {
            None => Err(ParseError::NeedMore),
            Some(pos_end) => Ok(pos_end + 1),
        }
    })?;

    let mut newlines = false;
    let seq_len = rb.mark_field(|buf: &[u8], eof: bool| {
        match memchr_both(b'\n', b'+', &buf) {
            (None, _) => Err(ParseError::NeedMore),
            (Some(pos_end), has_newline) => {
                newlines = has_newline;
                Ok(pos_end + 1)
            },
        }
    })?;

    let _ = rb.mark_field(|buf: &[u8], eof: bool| {
        // the quality header
        match memchr(b'\n', &buf) {
            None => Err(ParseError::NeedMore),
            Some(pos_end) => Ok(pos_end + 1),
        }
    })?;
        
    let _ = rb.mark_field(|buf: &[u8], eof: bool| {
        // the actual quality score
        if eof {
            Ok(min(seq_len, buf.len()))
        } else if seq_len > buf.len() {
            Err(ParseError::NeedMore)
        } else {
            Ok(seq_len)
        }
    })?;

    let mut fields: [&[u8]; 4];
    unsafe {
        fields = mem::uninitialized();
    }
    rb.fields(&mut fields);

    let id = match str::from_utf8(&fields[0][1..&fields[0].len() - 1]) { 
        Ok(v) => Ok(v),
        Err(_) => Err(ParseError::Invalid(String::from("FASTQ header not UTF8"))),
    }?;
    let seq = strip_whitespace(fields[1], newlines);
    let mut qual = fields[3];
    // remove newlines from the end of the quality record
    if qual[qual.len() - 1] == b'\n' {
        qual = &qual[..qual.len() - 1];
    }
    if qual[qual.len() - 1] == b'\r' {
        qual = &qual[..qual.len() - 1];
    }

    Ok((id, seq, Some(qual)))
}

// #[test]
// fn test_parsers() {
//     let parsed = fasta_record(b">\n", true);
//     let res = ("", b"".to_vec(), None);
//     assert_eq!(parsed, Ok((2, res)));
// 
//     let parsed = fasta_record(b">\n\n>", false);
//     let res = ("", b"".to_vec(), None);
//     assert_eq!(parsed, Ok((3, res)));
// 
//     let parsed = fasta_record(b">test\nagct\n>", false);
//     let res = ("test", b"agct".to_vec(), None);
//     assert_eq!(parsed, Ok((11, res)));
// 
//     let parsed = fasta_record(b">test2\nagct", true);
//     let res = ("test2", b"agct".to_vec(), None);
//     assert_eq!(parsed, Ok((11, res)));
// 
//     let parsed = fastq_record(b"@test\nagct\n+test\nAAAA\n", true);
//     let res = ("test", b"agct".to_vec(), Some(&b"AAAA"[..]));
//     assert_eq!(parsed, Ok((22, res)));
// }

/// Internal function abstracting over byte and file FASTX parsing
fn fastx_reader<'b, F, T>(reader: &'b mut T, first: u8, ref mut callback: F) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
          T: Read,
{
    let mut producer = RecBuffer::new(reader, 10000000usize);
    let parser = match first {
        b'>' => Ok(fasta_record as for<'a> fn(&'a mut RecBuffer) -> Result<SeqRecord<'a>, ParseError>),
        b'@' => Ok(fastq_record as for<'a> fn(&'a mut RecBuffer) -> Result<SeqRecord<'a>, ParseError>),
        _ => Err(ParseError::Invalid(String::from("Bad starting byte"))),
    }?;

    // TODO: replace this with some kind of futures event loop?
    loop {
        let record = parser(&mut producer);
        match record {
            Err(ParseError::EOF) => break,
            Ok(v) => callback(v),
            Err(e) => return Err(e),
        }
    }
    Ok(())
}

/// Parse a collection of bytes into FASTX records and calls `callback` on each.
///
pub fn fastx_bytes<'b, F>(bytes: &'b [u8], ref mut callback: F) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
{
    let mut cursor = Cursor::new(bytes);
    fastx_reader(&mut cursor, bytes[0], callback)
}


#[cfg(feature = "gz")]
pub fn fastx_file<F>(filename: &str, ref mut callback: F) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
{
    //! Parse a file (given its name) into FASTX records and calls `callback` on each.
    //!
    //! This version handles files compressed with gzip and should transparently
    //! decompress them.
    let mut f = File::open(&Path::new(filename))?;

    let mut first = vec![0];
    f.read(&mut first)?;
    if first[0] == 0x1F {
        let _ = f.seek(SeekFrom::Start(0));
        let mut gz_reader = GzDecoder::new(f)?;
        gz_reader.read(&mut first)?;
        // unfortunately, we can't seek back so we have to reopen everything
        drop(gz_reader);
        f = File::open(&Path::new(filename))?;
        gz_reader = GzDecoder::new(f)?;
        fastx_reader(&mut gz_reader, first[0], callback)
    } else {
        let _ = f.seek(SeekFrom::Start(0));
        fastx_reader(&mut f, first[0], callback)
    }
}

#[cfg(not(feature = "gz"))]
pub fn fastx_file<F>(filename: &str, ref mut callback: F) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
{
    //! Parse a file (given its name) into FASTX records and calls `callback` on each.
    let mut f = File::open(&Path::new(filename))?;

    let mut first = vec![0];
    f.read(&mut first)?;
    let _ = f.seek(SeekFrom::Start(0));
    fastx_reader(&mut f, first[0], callback)
}


#[test]
fn test_callback() {
    let mut i = 0;
    fastx_bytes(&b">test\nAGCT\n>test2\nGATC"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.0, "test");
                assert_eq!(&seq.1[..], &b"AGCT"[..]);
                assert_eq!(seq.2, None);
            },
            1 => {
                assert_eq!(seq.0, "test2");
                assert_eq!(&seq.1[..], &b"GATC"[..]);
                assert_eq!(seq.2, None);
            },
            _ => {
                assert!(false);
            },
        }
        i += 1;
    }).unwrap();
    assert_eq!(i, 2);

    i = 0;
    fastx_file("./tests/data/test.fa", |seq| {
        match i {
            0 => {
                assert_eq!(seq.0, "test");
                assert_eq!(&seq.1[..], b"AGCTGATCGA");
                assert_eq!(seq.2, None);
            },
            1 => {
                assert_eq!(seq.0, "test2");
                assert_eq!(&seq.1[..], b"TAGC");
                assert_eq!(seq.2, None);
            },
            _ => {
                assert!(false);
            },
        }
        i += 1;
    }).unwrap();
    assert_eq!(i, 2);
}

/// Given a SeqRecord and a quality cutoff, mask out low-quality bases with
/// `N` characters.
///
/// Experimental. Not currently used elsewhere in the codebase.
fn quality_mask<'a>(seq_rec: SeqRecord<'a>, ref score: u8) -> Vec<u8> {
    match seq_rec.2 {
        None => seq_rec.1.into_owned(),
        Some(quality) => {
            // could maybe speed this up by doing a copy of base and then
            // iterating though qual and masking?
            seq_rec.1
                   .iter()
                   .zip(quality.iter())
                   .map(|(base, qual)| {
                       if qual < score {
                           b'N'
                       } else {
                           base.clone()
                       }
                   })
                   .collect()

        },
    }
}

#[test]
fn test_quality_mask() {
    let seq_rec = ("", Cow::Borrowed(&b"AGCT"[..]), Some(&b"AAA0"[..]));
    let filtered = quality_mask(seq_rec, '5' as u8);
    assert_eq!(&filtered[..], &b"AGCN"[..]);
}


#[test]
fn test_fastq() {
    let mut i = 0;
    let res = fastx_bytes(&b"@test\nAGCT\n+test\n~~a!\n@test2\nTGCA\n+test\nWUI9"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.0, "test");
                assert_eq!(&seq.1[..], &b"AGCT"[..]);
                assert_eq!(seq.2, Some(&b"~~a!"[..]));
            },
            1 => {
                assert_eq!(seq.0, "test2");
                assert_eq!(&seq.1[..], &b"TGCA"[..]);
                assert_eq!(seq.2, Some(&b"WUI9"[..]));
            },
            _ => {
                assert!(false);
            },
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));
}


#[test]
fn test_wrapped_fasta() {
    let mut i = 0;
    let res = fastx_bytes(&b">test\nAGCT\nTCG\n>test2\nG"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.0, "test");
                assert_eq!(&seq.1[..], &b"AGCTTCG"[..]);
                assert_eq!(seq.2, None);
            },
            1 => {
                assert_eq!(seq.0, "test2");
                assert_eq!(&seq.1[..], &b"G"[..]);
                assert_eq!(seq.2, None);
            },
            _ => {
                assert!(false);
            },
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));
}
