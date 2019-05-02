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
use std::io::{stdin, Cursor, Read, Seek, SeekFrom};
use std::path::Path;
use std::str;

use memchr::memchr;

use crate::buffer::{FindRecord, RecBuffer, RecReader};
use crate::seq::SeqRecord;
use crate::util::{memchr_both, strip_whitespace, ParseError};

#[cfg(feature = "compression")]
use bzip2::read::BzDecoder;
#[cfg(feature = "compression")]
use flate2::read::MultiGzDecoder;
#[cfg(feature = "compression")]
use xz2::read::XzDecoder;
#[cfg(feature = "compression")]
use zip::read::ZipArchive;

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
        Some(Ok(FASTA { id: id, seq: seq }))
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
            // sometimes there are extra returns at the end of a file so we shouldn't blow up
            if buf[0] == b'\r' || buf[0] == b'\n' {
                return None;
            } else {
                return Some(Err(ParseError::Invalid(String::from("Bad FASTQ record"))));
            }
        }

        let id_end;
        match memchr(b'\n', &buf) {
            Some(i) => id_end = i + 1,
            None => return None,
        };
        let mut raw_id = &buf[1..id_end - 1];

        let seq_end;
        match memchr_both(b'\n', b'+', &buf[id_end..]) {
            Some(i) => seq_end = id_end + i + 1,
            None => return None,
        };
        let mut seq = &buf[id_end..seq_end - 1];

        let id2_end;
        match memchr(b'\n', &buf[seq_end..]) {
            Some(i) => id2_end = seq_end + i + 1,
            None => return None,
        };

        // we know the qual scores must be the same length as the sequence
        // so we can just do some arithmatic instead of memchr'ing
        let mut qual_end = id2_end + seq.len() + 1;
        let mut buffer_used = qual_end;
        if qual_end > buf.len() {
            if !self.last {
                // we need to pull more into the buffer
                return None;
            }
            // now do some math to figure out if the file doesn't end with a newline
            let windows_ending = if seq.last() == Some(&b'\r') { 1 } else { 0 };
            if !(qual_end == buf.len() + 1 + windows_ending) {
                return None;
            }
            buffer_used -= 1 + windows_ending;
            qual_end -= windows_ending;
        }
        let mut qual = &buf[id2_end..qual_end - 1];

        // clean up any extra '\r' from the id and seq
        if raw_id.len() > 0 && raw_id[raw_id.len() - 1] == b'\r' {
            raw_id = &raw_id[..raw_id.len() - 1];
            seq = &seq[..seq.len() - 1];
        }
        // we do qual separately in case this is the end of the file
        if qual.len() > 0 && qual[qual.len() - 1] == b'\r' {
            qual = &qual[..qual.len() - 1];
        }

        let id;
        match str::from_utf8(raw_id) {
            Ok(i) => id = i,
            Err(e) => return Some(Err(ParseError::from(e))),
        }
        self.pos += buffer_used;
        Some(Ok(FASTQ {
            id: id,
            seq: seq,
            qual: qual,
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
fn fastx_reader<F, R, T>(
    reader: &mut R,
    first_byte: Option<u8>,
    ref mut callback: F,
    type_callback: Option<&mut T>,
) -> Result<(), ParseError>
where
    F: for<'a> FnMut(SeqRecord<'a>) -> (),
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
where
    F: for<'a> FnMut(SeqRecord<'a>) -> (),
{
    let mut cursor = Cursor::new(bytes);
    fastx_reader(
        &mut cursor,
        None,
        callback,
        None::<&mut FnMut(&'static str) -> ()>,
    )
}

#[cfg(feature = "compression")]
pub fn fastx_stream<F, R, T>(
    mut reader: R,
    ref mut type_callback: T,
    ref mut callback: F,
) -> Result<(), ParseError>
where
    F: for<'a> FnMut(SeqRecord<'a>) -> (),
    R: Read + Seek,
    T: FnMut(&'static str) -> (),
{
    //! Opens a `Read` stream and parses the FASTX records out. Also takes a "type_callback"
    //! that gets called as soon as we determine if the records are FASTA or FASTQ.
    //!  If a file starts with a gzip or other header, transparently decompress it.
    let mut first = vec![0];
    reader.read(&mut first)?;
    if first[0] == 0x1F {
        // gz files
        reader.read(&mut first)?;
        if first[0] != 0x8B {
            return Err(ParseError::Invalid(String::from("Bad starting bytes")));
        }
        let _ = reader.seek(SeekFrom::Start(0));
        let mut gz_reader = MultiGzDecoder::new(reader);
        fastx_reader(&mut gz_reader, None, callback, Some(type_callback))
    } else if first[0] == 0x42 {
        // bz files
        reader.read(&mut first)?;
        if first[0] != 0x5A {
            return Err(ParseError::Invalid(String::from("Bad starting bytes")));
        }
        let _ = reader.seek(SeekFrom::Start(0));
        let mut bz_reader = BzDecoder::new(reader);
        fastx_reader(&mut bz_reader, None, callback, Some(type_callback))
    } else if first[0] == 0xFD {
        // xz files
        reader.read(&mut first)?;
        if first[0] != 0x37 {
            return Err(ParseError::Invalid(String::from("Bad starting bytes")));
        }
        let _ = reader.seek(SeekFrom::Start(0));
        let mut xz_reader = XzDecoder::new(reader);
        fastx_reader(&mut xz_reader, None, callback, Some(type_callback))
    } else if first[0] == 0x50 {
        // zip files
        reader.read(&mut first)?;
        if first[0] != 0x4b {
            return Err(ParseError::Invalid(String::from("Bad starting bytes")));
        }
        let _ = reader.seek(SeekFrom::Start(0));

        let mut zip_archive = ZipArchive::new(reader)?;
        if zip_archive.len() != 1 {
            return Err(ParseError::Invalid(String::from(
                "Zip archives has more than one file",
            )));
        }
        let mut zip_reader = zip_archive.by_index(0)?;
        fastx_reader(&mut zip_reader, None, callback, Some(type_callback))
    } else {
        fastx_reader(&mut reader, Some(first[0]), callback, Some(type_callback))
    }
}

#[cfg(feature = "compression")]
pub fn fastx_cli<F, T>(
    filename: &str,
    ref mut type_callback: T,
    ref mut callback: F,
) -> Result<(), ParseError>
where
    F: for<'a> FnMut(SeqRecord<'a>) -> (),
    T: FnMut(&'static str) -> (),
{
    //! Opens files (or stdin, if a dash is provided instead) and reads FASTX records out. Also
    //! takes a "type_callback" that gets called as soon as we determine if the records are FASTA
    //! or FASTQ.  If a file starts with a gzip or other header, transparently decompress it.
    if filename == "-" {
        let sin = stdin();
        let mut lock = sin.lock();
        return fastx_reader(&mut lock, None, callback, Some(type_callback));
    }

    let f = File::open(&Path::new(filename))?;
    fastx_stream(f, type_callback, callback)
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
            _ => unreachable!("Too many records"),
        }
        i += 1;
    });
    assert_eq!(res, Ok(()));
    assert_eq!(i, 2);

    i = 0;
    let res = fastx_cli(
        "./tests/data/test.fa",
        |filetype| {
            assert_eq!(filetype, "FASTA");
        },
        |seq| {
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
                _ => unreachable!("Too many records"),
            }
            i += 1;
        },
    );
    assert_eq!(res, Ok(()));
    assert_eq!(i, 2);

    let res = fastx_cli(
        "./tests/data/bad_test.fa",
        |_| {
            unreachable!("This is not a valid file type");
        },
        |_| {
            unreachable!("No valid records in this file to parse");
        },
    );
    assert_eq!(
        res,
        Err(ParseError::Invalid(String::from("Bad starting byte")))
    );
}

#[cfg(feature = "compression")]
#[test]
fn test_compressed() {
    let test_files = [
        "./tests/data/test.fa.gz",
        "./tests/data/test.fa.bz2",
        "./tests/data/test.fa.xz",
        "./tests/data/test.fa.zip",
    ];

    for test_file in test_files.iter() {
        let mut i = 0;
        let res = fastx_cli(
            test_file,
            |filetype| {
                assert_eq!(filetype, "FASTA");
            },
            |seq| {
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
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(res, Ok(()));
        assert_eq!(i, 2);
    }
}

#[test]
fn test_fastq() {
    let mut i = 0;
    let res = fastx_bytes(
        &b"@test\nAGCT\n+test\n~~a!\n@test2\nTGCA\n+test\nWUI9"[..],
        |seq| {
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
                _ => unreachable!("Too many records"),
            }
            i += 1;
        },
    );
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));

    let mut i = 0;
    let res = fastx_bytes(
        &b"@test\r\nAGCT\r\n+test\r\n~~a!\r\n@test2\r\nTGCA\r\n+test\r\nWUI9"[..],
        |seq| {
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
                _ => unreachable!("Too many records"),
            }
            i += 1;
        },
    );
    assert_eq!(res, Ok(()));
    assert_eq!(i, 2);
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
            _ => unreachable!("Too many records"),
        }
        i += 1;
    });
    assert_eq!(res, Ok(()));
    assert_eq!(i, 2);

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
            _ => unreachable!("Too many records"),
        }
        i += 1;
    });
    assert_eq!(res, Ok(()));
    assert_eq!(i, 2);
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
            _ => unreachable!("Too many records"),
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
            _ => unreachable!("Too many records"),
        }
        i += 1;
    });
    assert_eq!(i, 1);
    assert_eq!(res, Err(ParseError::PrematureEOF));

    // we allow a few extra newlines at the ends of FASTQs
    let mut i = 0;
    let res = fastx_bytes(&b"@test\nAGCT\n+test\n~~a!\n\n"[..], |seq| {
        match i {
            0 => {
                assert_eq!(seq.id, "test");
                assert_eq!(&seq.seq[..], &b"AGCT"[..]);
                assert_eq!(&seq.qual.unwrap()[..], &b"~~a!"[..]);
            },
            _ => unreachable!("Too many records"),
        }
        i += 1;
    });
    assert_eq!(i, 1);
    assert_eq!(res, Ok(()));

    // but if there's additional data past the newlines it's an error
    let mut i = 0;
    let res = fastx_bytes(
        &b"@test\nAGCT\n+test\n~~a!\n\n@TEST\nA\n+TEST\n~"[..],
        |seq| {
            match i {
                0 => {
                    assert_eq!(seq.id, "test");
                    assert_eq!(&seq.seq[..], &b"AGCT"[..]);
                    assert_eq!(&seq.qual.unwrap()[..], &b"~~a!"[..]);
                },
                _ => unreachable!("Too many records"),
            }
            i += 1;
        },
    );
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
            _ => unreachable!("Too many records"),
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
            _ => unreachable!("Too many records"),
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
            _ => unreachable!("Too many records"),
        }
        i += 1;
    });
    assert_eq!(i, 2);
    assert_eq!(res, Ok(()));
}

#[test]
fn test_buffer() {
    let mut buf: RecBuffer<FASTA> = RecBuffer::from_bytes(b">test\nACGT");
    let rec = buf.next().unwrap().unwrap();
    assert_eq!(rec.id, "test", "Record has the right ID");
    assert_eq!(rec.seq, b"ACGT", "Record has the right sequence");

    let mut buf: RecBuffer<FASTA> = RecBuffer::from_bytes(b">test");
    assert!(buf.next().is_none(), "Incomplete record returns None");
}

#[test]
fn test_fastq_across_buffer() {
    let test_seq = b"@A\nA\n+A\nA\n@B\nA\n+B\n!";
    let mut cursor = Cursor::new(test_seq);
    // the buffer is aligned to the first record
    let mut rec_reader = RecReader::new(&mut cursor, 9, b"").unwrap();

    let used = {
        let mut rec_buffer = rec_reader.get_buffer::<FASTQ>();
        for s in rec_buffer.by_ref() {
            // record is incomplete
            panic!("No initial record should be parsed")
        }
        rec_buffer.pos
    };

    // refill the buffer, but we're not done quite yet
    assert_eq!(rec_reader.refill(used).unwrap(), false);

    // now we should see both records
    let mut rec_buffer = rec_reader.get_buffer::<FASTQ>();

    // there should be a record assuming the parser
    // handled the buffer boundary
    let iterated_seq = rec_buffer.by_ref().next();
    let seq = iterated_seq.unwrap();
    assert_eq!(seq.unwrap().id, "A");

    // but not another because the buffer's too short
    let iterated_seq = rec_buffer.by_ref().next();
    assert!(iterated_seq.is_none());

    // TODO: refill and check for the last record
}
