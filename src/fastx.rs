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
use std::fs::File;
use std::io::{Cursor, Read, Seek, SeekFrom, stdin};
use std::path::Path;
use std::str;

use memchr::memchr;

use buffer::{parse, RecBuffer};
use seq::SeqRecord;
use util::{ParseError, memchr_both};

#[cfg(feature = "gz")]
use flate2::read::GzDecoder;


// TODO: make `id` fields &str
#[derive(Debug)]
struct FASTA<'a> {
    id: &'a [u8],
    seq: &'a [u8],
}


#[derive(Debug)]
struct FASTQ<'a> {
    id: &'a [u8],
    seq: &'a [u8],
    qual: &'a [u8],
}


impl<'a> Iterator for RecBuffer<'a, FASTA<'static>> {
    type Item = Result<FASTA<'a>, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let start = self.pos;
        let id_end;
        match memchr(b'\n', &self.buf[start..]) {
            Some(i) => id_end = i,
            None => return None,
        };

        let seq_end;
        match (memchr_both(b'\n', b'>', &self.buf[id_end..]).0, self.last) {
            (Some(i), _) => seq_end = i,
            (None, true) => seq_end = self.buf.len(),
            (None, false) => return None,
        };
        if self.buf[start] != b'>' {
            return Some(Err(ParseError::Invalid(String::from("Bad FASTA record"))));
        }
        self.pos = seq_end;
        Some(Ok(FASTA {
            id: &self.buf[start..id_end],
            seq: &self.buf[id_end..seq_end],
        }))
    }
}


impl<'a> Iterator for RecBuffer<'a, FASTQ<'a>> {
    type Item = Result<FASTQ<'a>, ParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let start = self.pos;

        let id_end;
        match memchr(b'\n', &self.buf[start..]) {
            Some(i) => id_end = i,
            None => return None,
        };
        let seq_end;
        match memchr_both(b'\n', b'+', &self.buf[id_end..]).0 {
            Some(i) => seq_end = i,
            None => return None,
        };
        let id2_end;
        match memchr(b'\n', &self.buf[seq_end..]) {
            Some(i) => id2_end = i,
            None => return None,
        };
        // we know the qual scores must be the same length as the sequence
        // so we can just do some arithmatic instead of memchr'ing
        if (self.buf.len() - id2_end) < (id_end - start) {
            return None;
        }
        let qual_end = id2_end + id_end - start;
        if self.buf[start] != b'@' {
            return Some(Err(ParseError::Invalid(String::from("Bad FASTA record"))));
        }
        self.pos = qual_end;
        Some(Ok(FASTQ {
            id: &self.buf[start..id_end],
            seq: &self.buf[id_end..seq_end],
            qual: &self.buf[id2_end..qual_end],
        }))
    }
}

impl<'a> From<FASTA<'a>> for SeqRecord<'a> {
    fn from(fasta: FASTA<'a>) -> SeqRecord<'a> {
        SeqRecord::new(str::from_utf8(fasta.id).unwrap(), Cow::from(fasta.seq), None)
    }

}

impl<'a> From<FASTQ<'a>> for SeqRecord<'a> {
    fn from(fastq: FASTQ<'a>) -> SeqRecord<'a> {
        SeqRecord::new(str::from_utf8(fastq.id).unwrap(), Cow::from(fastq.seq), Some(fastq.qual))
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
    match first[0] {
        b'>' => {
            if let Some(f) = type_callback {
                f("FASTA");
            }
            parse::<FASTA, ParseError, _>(reader, &[b'>'], |record| {
                callback(SeqRecord::from(record));
                Ok(())
            });
        },
        b'@' => {
            if let Some(f) = type_callback {
                f("FASTQ");
            }
            // parse::<FASTQ>(reader, &b'@', |record| {
            //     callback(record.to_seq()?);
            //     Ok(())
            // });
        },
        _ => return Err(ParseError::Invalid(String::from("Bad starting byte"))),
    };
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
        let mut gz_reader = GzDecoder::new(reader)?;
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
        let mut gz_reader = GzDecoder::new(f)?;
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
        let mut gz_reader = GzDecoder::new(f)?;
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
                assert_eq!(seq.qual, Some(&b"~~a!"[..]));
            },
            1 => {
                assert_eq!(seq.id, "test2");
                assert_eq!(&seq.seq[..], &b"TGCA"[..]);
                assert_eq!(seq.qual, Some(&b"WUI9"[..]));
            },
            _ => assert!(false),
        }
        i += 1;
    });
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
    assert_eq!(res, Err(ParseError::PrematureEOF));

    let mut i = 0;
    let res = fastx_bytes(&b"@test\nAGCT\n+test\n~~a!\n@test2\nTGCA"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], &b"AGCT"[..]);
                assert_eq!(seq.qual, Some(&b"~~a!"[..]));
            },
            _ => assert!(false),
        }
        i += 1;
    });
    assert_eq!(res, Err(ParseError::PrematureEOF));
}
