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
use std::io::{Cursor, Read, Seek, SeekFrom, stdin};
use std::path::Path;
use std::str;

use memchr::memchr;

use buffer::{RecBuffer, RecReader, FindRecord};
use seq::SeqRecord;
use util::{ParseError, memchr_both, strip_whitespace};

#[cfg(feature = "gz")]
use flate2::read::GzDecoder;


#[derive(Debug)]
struct FASTA<'a> {
    id: &'a str,
    seq: &'a [u8],
}


#[derive(Debug)]
struct FASTQ<'a> {
    id: &'a str,
    seq: &'a [u8],
    qual: &'a [u8],
}


impl<'a> Iterator for RecBuffer<'a, FASTA<'static>> {
    type Item = Result<FASTA<'a>, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let buf = &self.buf[self.pos..];
        if buf.len() == 0 {
            return None;
        }
        if buf[0] != b'>' {
            return Some(Err(ParseError::Invalid(String::from("Bad FASTA record"))));
        }

        let id_end;
        match memchr(b'\n', &buf) {
            Some(i) => id_end = i + 1,
            None => return None,
        };
        let mut raw_id = &buf[1..id_end - 1];
        if raw_id.len() > 0 && raw_id[raw_id.len() - 1] == b'\r' {
            raw_id = &raw_id[..raw_id.len() - 1];
        }
        let id;
        match str::from_utf8(raw_id) {
            Ok(i) => id = i,
            Err(e) => return Some(Err(ParseError::from(e))),
        }

        let seq_end;
        match (memchr_both(b'\n', b'>', &buf[id_end..]), self.last) {
            (Some(i), _) => seq_end = id_end + i + 1,
            (None, true) => seq_end = self.buf.len(),
            (None, false) => return None,
        };
        let mut seq = &buf[id_end..seq_end];
        if seq[seq.len() - 1] == b'\r' {
            seq = &seq[..seq.len()];
        }

        self.pos += seq_end;
        Some(Ok(FASTA {
            id: id,
            seq: seq,
        }))
    }
}


impl<'a> Iterator for RecBuffer<'a, FASTQ<'a>> {
    type Item = Result<FASTQ<'a>, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pos >= self.buf.len() {
            return None;
        }
        let buf = &self.buf[self.pos..];
        if buf[0] != b'@' {
            return Some(Err(ParseError::Invalid(String::from("Bad FASTQ record"))));
        }

        let id_end;
        match memchr(b'\n', &buf) {
            Some(i) => id_end = i + 1,
            None => return None,
        };

        let seq_end;
        match memchr_both(b'\n', b'+', &buf[id_end..]) {
            Some(i) => seq_end = id_end + i,
            None => return None,
        };
        let mut seq = &buf[id_end..seq_end];
        if seq.len() > 0 && seq[seq.len() - 1] == b'\r' {
            seq = &seq[..seq.len() - 1];
        }

        let id2_end;
        match memchr(b'\n', &buf[seq_end + 1..]) {
            Some(i) => id2_end = seq_end + i + 1,
            None => return None,
        };
        // we know the qual scores must be the same length as the sequence
        // so we can just do some arithmatic instead of memchr'ing
        if (buf.len() - id2_end) < seq.len() {
            return None;
        }
        let mut raw_id = &buf[1..id_end - 1];
        if raw_id.len() > 0 && raw_id[raw_id.len() - 1] == b'\r' {
            raw_id = &raw_id[..raw_id.len() - 1];
        }
        let id;
        match str::from_utf8(raw_id) {
            Ok(i) => id = i,
            Err(e) => return Some(Err(ParseError::from(e))),
        }
        // if there are \r or \n we need to skip them still
        // FIXME: 2 isn't always right
        // println!("{:?}", &buf[self.pos..self.pos+id2_end + seq_end - id_end + 2]);
        self.pos += min(id2_end + seq_end - id_end + 2, buf.len());
        Some(Ok(FASTQ {
            id: id,
            seq: seq,
            qual: &buf[id2_end + 1..id2_end + 1 + seq.len()],
        }))
    }
}


fn is_finished<T>(rb: &RecBuffer<T>) -> bool {
    if !rb.last {
        return false;
    }
    if rb.pos == rb.buf.len() {
        return true;
    }
    for c in &rb.buf[rb.pos..] {
        if c != &b'\r' && c != &b'\n' {
            return false;
        }
    }
    true
}


impl<'a> FindRecord for RecBuffer<'a, FASTA<'a>> {
    fn move_to_next(&mut self) {
        unimplemented!("");
    }

    fn is_finished(&self) -> bool {
        is_finished(&self)
    }
}


impl<'a> FindRecord for RecBuffer<'a, FASTQ<'a>> {
    fn move_to_next(&mut self) {
        unimplemented!("");
    }

    fn is_finished(&self) -> bool {
        is_finished(&self)
    }
}



impl<'a> From<FASTA<'a>> for SeqRecord<'a> {
    fn from(fasta: FASTA<'a>) -> SeqRecord<'a> {
        SeqRecord::new(fasta.id, Cow::from(strip_whitespace(fasta.seq)), None)
    }
}

impl<'a> From<FASTQ<'a>> for SeqRecord<'a> {
    fn from(fastq: FASTQ<'a>) -> SeqRecord<'a> {
        SeqRecord::new(fastq.id, Cow::from(fastq.seq), Some(fastq.qual))
    }
}


/// Internal function abstracting over byte and file FASTX parsing
fn fastx_reader<F, R, T>(reader: &mut R, first_byte: Option<u8>, ref mut callback: F, type_callback: Option<&mut T>) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
          R: Read,
          T: ?Sized + FnMut(&'static str) -> (),
{
    let mut first = vec![0];
    match first_byte {
        Some(b) => first[0] = b,
        None => {
            reader.read(&mut first)?;
        },
    }
    let mut rec_reader = RecReader::new(reader, 10_000_000, &first)?;
    match first[0] {
        b'>' => {
            if let Some(f) = type_callback {
                f("FASTA");
            }
            loop {
                let used = {
                    let mut rec_buffer = rec_reader.get_buffer::<FASTA>();
                    for s in rec_buffer.by_ref() {
                        callback(SeqRecord::from(s?));
                    }
                    rec_buffer.pos
                };
                if rec_reader.refill(used)? {
                    break;
                }
            }
        },
        b'@' => {
            if let Some(f) = type_callback {
                f("FASTQ");
            }
            loop {
                let used = {
                    let mut rec_buffer = rec_reader.get_buffer::<FASTQ>();
                    for s in rec_buffer.by_ref() {
                        callback(SeqRecord::from(s?));
                    }
                    rec_buffer.pos
                };
                if rec_reader.refill(used)? {
                    break;
                }
            }
        },
        _ => return Err(ParseError::Invalid(String::from("Bad starting byte"))),
    };
    // check if there's anything left stuff in the buffer (besides returns)
    let rec_buffer = rec_reader.get_buffer::<FASTA>();
    if !rec_buffer.is_finished() {
        return Err(ParseError::PrematureEOF);
    }
    Ok(())
}


/// Parse a array of bytes into FASTX records and calls `callback` on each.
pub fn fastx_bytes<'b, F>(bytes: &'b [u8], ref mut callback: F) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
{
    let mut cursor = Cursor::new(bytes);
    fastx_reader(&mut cursor, None, callback, None::<&mut FnMut(&'static str) -> ()>)
}


/// Parse stdin into FASTX records and call `callback` on each.
#[deprecated(since="0.1.3", note="please use `fastx_stream` instead")]
pub fn fastx_stdin<F>(ref mut callback: F) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> ()
{
    let sin = stdin();
    let mut lock = sin.lock();
    fastx_reader(&mut lock, None, callback, None::<&mut FnMut(&'static str) -> ()>)
}


#[cfg(feature = "gz")]
pub fn fastx_stream<F, R, T>(mut reader: R, ref mut type_callback: T, ref mut callback: F) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
          R: Read + Seek,
          T: FnMut(&'static str) -> (),
{
    //! Opens a `Read` stream and parses the FASTX records out. Also takes a "type_callback"
    //! that gets called as soon as we determine if the records are FASTA or FASTQ.
    //!  If a file starts with a gzip header, transparently decompress it.
    let mut first = vec![0];
    reader.read(&mut first)?;
    if first[0] == 0x1F {
        reader.read(&mut first)?;
        if first[0] != 0x8B {
            return Err(ParseError::Invalid(String::from("Bad starting bytes")));
        }
        let _ = reader.seek(SeekFrom::Start(0));
        let mut gz_reader = GzDecoder::new(reader);
        fastx_reader(&mut gz_reader, None, callback, Some(type_callback))
    } else {
        fastx_reader(&mut reader, Some(first[0]), callback, Some(type_callback))
    }
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
        f.read(&mut first)?;
        if first[0] != 0x8B {
            return Err(ParseError::Invalid(String::from("Bad starting bytes")));
        }
        let _ = f.seek(SeekFrom::Start(0));
        let mut gz_reader = GzDecoder::new(f);
        fastx_reader(&mut gz_reader, None, callback, None::<&mut FnMut(&'static str) -> ()>)
    } else {
        fastx_reader(&mut f, Some(first[0]), callback, None::<&mut FnMut(&'static str) -> ()>)
    }
}

#[cfg(not(feature = "gz"))]
pub fn fastx_file<F>(filename: &str, ref mut callback: F) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
{
    //! Parse a file (given its name) into FASTX records and calls `callback` on each.
    let mut f = File::open(&Path::new(filename))?;

    fastx_reader(&mut f, None, callback, None::<&mut FnMut(&'static str) -> ()>)
}


#[cfg(feature = "gz")]
pub fn fastx_cli<F, T>(filename: &str, ref mut type_callback: T, ref mut callback: F) -> Result<(), ParseError>
    where F: for<'a> FnMut(SeqRecord<'a>) -> (),
          T: FnMut(&'static str) -> (),
{
    //! Opens files (or stdin, if a dash is provided instead) and reads FASTX records out. Also
    //! takes a "type_callback" that gets called as soon as we determine if the records are FASTA
    //! or FASTQ.  If a file starts with a gzip header, transparently decompress it.
    if filename == "-" {
        let sin = stdin();
        let mut lock = sin.lock();
        return fastx_reader(&mut lock, None, callback, Some(type_callback));
    }

    let mut f = File::open(&Path::new(filename))?;

    let mut first = vec![0];
    f.read(&mut first)?;
    if first[0] == 0x1F {
        f.read(&mut first)?;
        if first[0] != 0x8B {
            return Err(ParseError::Invalid(String::from("Bad starting bytes")));
        }
        let _ = f.seek(SeekFrom::Start(0));
        let mut gz_reader = GzDecoder::new(f);
        fastx_reader(&mut gz_reader, None, callback, Some(type_callback))
    } else {
        fastx_reader(&mut f, Some(first[0]), callback, Some(type_callback))
    }
}


#[test]
fn test_callback() {
    let mut i = 0;
    let res = fastx_bytes(&b">test\nAGCT\n>test2\nGATC"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], &b"AGCT"[..]);
                assert_eq!(seq.qual, None);
            },
            1 => {
                assert_eq!(seq.id, "test2");
                assert_eq!(&seq.seq[..], &b"GATC"[..]);
                assert_eq!(seq.qual, None);
            },
            _ => {
                assert!(false);
            },
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));

    i = 0;
    let res = fastx_file("./tests/data/test.fa", |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], b"AGCTGATCGA");
                assert_eq!(seq.qual, None);
            },
            1 => {
                assert_eq!(seq.id, "test2");
                assert_eq!(&seq.seq[..], b"TAGC");
                assert_eq!(seq.qual, None);
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));
}

#[test]
fn test_fastq() {
    let mut i = 0;
    let res = fastx_bytes(&b"@test\nAGCT\n+test\n~~a!\n@test2\nTGCA\n+test\nWUI9"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], &b"AGCT"[..]);
                assert_eq!(&seq.qual.unwrap()[..], &b"~~a!"[..]);
            },
            1 => {
                assert_eq!(seq.id, "test2");
                assert_eq!(&seq.seq[..], &b"TGCA"[..]);
                assert_eq!(&seq.qual.unwrap()[..], &b"WUI9"[..]);
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));

    let mut i = 0;
    let res = fastx_bytes(&b"@test\r\nAGCT\r\n+test\r\n~~a!\r\n@test2\r\nTGCA\r\n+test\r\nWUI9"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], &b"AGCT"[..]);
                assert_eq!(&seq.qual.unwrap()[..], &b"~~a!"[..]);
            },
            1 => {
                assert_eq!(seq.id, "test2");
                assert_eq!(&seq.seq[..], &b"TGCA"[..]);
                assert_eq!(&seq.qual.unwrap()[..], &b"WUI9"[..]);
            },
            _ => assert!(false),
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
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], &b"AGCTTCG"[..]);
                assert_eq!(seq.qual, None);
            },
            1 => {
                assert_eq!(seq.id, "test2");
                assert_eq!(&seq.seq[..], &b"G"[..]);
                assert_eq!(seq.qual, None);
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));

    let mut i = 0;
    let res = fastx_bytes(&b">test\r\nAGCT\r\nTCG\r\n>test2\r\nG"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], &b"AGCTTCG"[..]);
                assert_eq!(seq.qual, None);
            },
            1 => {
                assert_eq!(seq.id, "test2");
                assert_eq!(&seq.seq[..], &b"G"[..]);
                assert_eq!(seq.qual, None);
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));
}

#[test]
fn test_premature_endings() {
    let mut i = 0;
    let res = fastx_bytes(&b">test\nAGCT\n>test2"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], &b"AGCT"[..]);
                assert_eq!(seq.qual, None);
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(i, 1);
    assert_eq!(res, Err(ParseError::PrematureEOF));

    let mut i = 0;
    let res = fastx_bytes(&b"@test\nAGCT\n+test\n~~a!\n@test2\nTGCA"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], &b"AGCT"[..]);
                assert_eq!(&seq.qual.unwrap()[..], &b"~~a!"[..]);
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(i, 1);
    assert_eq!(res, Err(ParseError::PrematureEOF));
}

#[test]
fn test_empty_records() {
    let mut i = 0;
    let res = fastx_bytes(&b"@\n\n+\n\n@test2\nTGCA\n+test2\n~~~~\n"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "");
                assert_eq!(&seq.seq[..], &b""[..]);
                assert_eq!(&seq.qual.unwrap()[..], &b""[..]);
            },
            1 => {
                assert_eq!(seq.id, "test2");
                assert_eq!(&seq.seq[..], &b"TGCA"[..]);
                assert_eq!(&seq.qual.unwrap()[..], &b"~~~~"[..]);
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));

    let mut i = 0;
    let res = fastx_bytes(&b">\n\n>shine\nAGGAGGU"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "");
                assert_eq!(&seq.seq[..], &b""[..]);
                assert_eq!(seq.qual, None);
            },
            1 => {
                assert_eq!(seq.id, "shine");
                assert_eq!(&seq.seq[..], &b"AGGAGGU"[..]);
                assert_eq!(seq.qual, None);
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));

    let mut i = 0;
    let res = fastx_bytes(&b">\r\n\r\n>shine\r\nAGGAGGU"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "");
                assert_eq!(&seq.seq[..], &b""[..]);
                assert_eq!(seq.qual, None);
            },
            1 => {
                assert_eq!(seq.id, "shine");
                assert_eq!(&seq.seq[..], &b"AGGAGGU"[..]);
                assert_eq!(seq.qual, None);
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));
}
