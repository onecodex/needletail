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

mod fasta;
mod fastq;

use std::cmp::min;
use std::io::{Cursor, Read};
use std::str;

#[cfg(feature = "compression")]
use bzip2::read::BzDecoder;
#[cfg(feature = "compression")]
use flate2::read::MultiGzDecoder;
#[cfg(feature = "compression")]
use xz2::read::XzDecoder;

use crate::buffer::RecReader;
pub use crate::formats::fasta::FASTA;
pub use crate::formats::fastq::FASTQ;
use crate::seq::Sequence;
use crate::util::{ParseError, ParseErrorType};


/// Internal function abstracting over byte and file FASTX parsing
fn seq_reader<F, R, T>(
    reader: &mut R,
    mut callback: F,
    type_callback: &mut T,
) -> Result<(), ParseError>
where
    F: for<'a> FnMut(Sequence<'a>) -> (),
    R: Read,
    T: ?Sized + FnMut(&'static str) -> (),
{
    let mut first = vec![0];
    reader.read_exact(&mut first)?;
    match first[0] {
        b'>' => type_callback("FASTA"),
        b'@' => type_callback("FASTQ"),
        _ => (),
    }
    let mut rec_reader = RecReader::new(reader, 10_000_000, &first)?;
    let mut record_count = 0;
    loop {
        let used = match first[0] {
            b'>' => {
                let mut rec_buffer = rec_reader.get_buffer::<FASTA>(record_count);
                for s in rec_buffer.by_ref() {
                    callback(Sequence::from(s?));
                }
                record_count += rec_buffer.count;
                rec_buffer.pos
            },
            b'@' => {
                let mut rec_buffer = rec_reader.get_buffer::<FASTQ>(record_count);
                for s in rec_buffer.by_ref() {
                    callback(Sequence::from(s?));
                }
                record_count += rec_buffer.count;
                rec_buffer.pos
            },
            _ => {
                return Err(
                    ParseError::new("Bad starting byte", ParseErrorType::InvalidHeader)
                        .record(0)
                        .context(String::from_utf8_lossy(&first)),
                )
            },
        };
        if rec_reader.refill(used)? {
            break;
        }
    }
    // check if there's anything left stuff in the buffer (besides returns)
    let rec_buffer = rec_reader.get_buffer::<FASTA>(record_count);
    if !rec_buffer.last {
        return Err(
            ParseError::new("File ended abruptly", ParseErrorType::PrematureEOF)
                .record(record_count),
        );
    }
    for c in &rec_buffer.buf[rec_buffer.pos..] {
        if c != &b'\r' && c != &b'\n' {
            let end = min(rec_buffer.pos + 16, rec_buffer.buf.len());
            let context = String::from_utf8_lossy(&rec_buffer.buf[rec_buffer.pos..end]);
            return Err(ParseError::new(
                "File had extra data past end of records",
                ParseErrorType::PrematureEOF,
            )
            .record(record_count)
            .context(context));
        }
    }
    Ok(())
}

#[cfg(not(feature = "compression"))]
pub fn parse_sequences<F, R, T>(
    mut reader: R,
    mut type_callback: T,
    callback: F,
) -> Result<(), ParseError>
where
    F: for<'a> FnMut(Sequence<'a>) -> (),
    R: Read,
    T: FnMut(&'static str) -> (),
{
    //! Opens a `Read` stream and parses the FASTX records out. Also takes a "type_callback"
    //! that gets called as soon as we determine if the records are FASTA or FASTQ.
    fastx_reader(&mut reader, callback, &mut type_callback)
}

#[cfg(feature = "compression")]
pub fn parse_sequences<F, R, T>(
    mut reader: R,
    mut type_callback: T,
    callback: F,
) -> Result<(), ParseError>
where
    F: for<'a> FnMut(Sequence<'a>) -> (),
    R: Read,
    T: FnMut(&'static str) -> (),
{
    //! Opens a `Read` stream and parses the FASTX records out. Also takes a "type_callback"
    //! that gets called as soon as we determine if the records are FASTA or FASTQ.
    //!  If a file starts with a gzip or other header, transparently decompress it.
    let mut first = vec![0];
    reader.read_exact(&mut first)?;
    if first[0] == 0x1F {
        // gz files
        reader.read_exact(&mut first)?;
        if first[0] != 0x8B {
            return Err(ParseError::new(
                "Bad gz header",
                ParseErrorType::BadCompression,
            ));
        }
        let cursor = Cursor::new(vec![0x1F, 0x8B]);
        let mut gz_reader = MultiGzDecoder::new(cursor.chain(reader));

        seq_reader(&mut gz_reader, callback, &mut type_callback)
    } else if first[0] == 0x42 {
        // bz files
        reader.read_exact(&mut first)?;
        if first[0] != 0x5A {
            return Err(ParseError::new(
                "Bad bz header",
                ParseErrorType::BadCompression,
            ));
        }
        let cursor = Cursor::new(vec![0x42, 0x5A]);
        let mut bz_reader = BzDecoder::new(cursor.chain(reader));

        seq_reader(&mut bz_reader, callback, &mut type_callback)
    } else if first[0] == 0xFD {
        // xz files
        reader.read_exact(&mut first)?;
        if first[0] != 0x37 {
            return Err(ParseError::new(
                "Bad xz header",
                ParseErrorType::BadCompression,
            ));
        }
        let cursor = Cursor::new(vec![0xFD, 0x37]);
        let mut xz_reader = XzDecoder::new(cursor.chain(reader));

        seq_reader(&mut xz_reader, callback, &mut type_callback)
    } else {
        let cursor = Cursor::new(first);
        let mut reader = cursor.chain(reader);
        seq_reader(&mut reader, callback, &mut type_callback)
    }
}

#[cfg(test)]
mod test {
    use std::fs::File;
    use std::io::Cursor;
    use std::path::Path;

    use crate::buffer::{RecBuffer, RecReader};
    use crate::formats::{parse_sequences, FASTA, FASTQ};
    use crate::util::ParseErrorType;

    fn seq(s: &[u8]) -> Cursor<&[u8]> {
        Cursor::new(&s[..])
    }

    #[test]
    fn test_callback() {
        let mut i = 0;
        let res = parse_sequences(
            seq(b">test\nAGCT\n>test2\nGATC"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCT");
                        assert_eq!(seq.qual, None);
                    },
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"GATC");
                        assert_eq!(seq.qual, None);
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(res, Ok(()));
        assert_eq!(i, 2);

        i = 0;
        let file = File::open(&Path::new("./tests/data/test.fa")).unwrap();
        let res = parse_sequences(
            file,
            |filetype| {
                assert_eq!(filetype, "FASTA");
            },
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCTGATCGA");
                        assert_eq!(seq.qual, None);
                    },
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
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

        let file = File::open(&Path::new("./tests/data/bad_test.fa")).unwrap();
        let res = parse_sequences(
            file,
            |_| {
                unreachable!("This is not a valid file type");
            },
            |_| {
                unreachable!("No valid records in this file to parse");
            },
        );
        let e = res.unwrap_err();
        assert_eq!(e.error_type, ParseErrorType::InvalidHeader);
        assert_eq!(e.record, 0);
        assert_eq!(e.msg, String::from("Bad starting byte"));
    }

    #[cfg(feature = "compression")]
    #[test]
    fn test_compressed() {
        let test_files = [
            "./tests/data/test.fa.gz",
            "./tests/data/test.fa.bz2",
            "./tests/data/test.fa.xz",
        ];

        for test_file in test_files.iter() {
            let mut i = 0;
            let file = File::open(&Path::new(test_file)).unwrap();
            let res = parse_sequences(
                file,
                |filetype| {
                    assert_eq!(filetype, "FASTA");
                },
                |seq| {
                    match i {
                        0 => {
                            assert_eq!(&seq.id[..], b"test");
                            assert_eq!(&seq.seq[..], b"AGCTGATCGA");
                            assert_eq!(seq.qual, None);
                        },
                        1 => {
                            assert_eq!(&seq.id[..], b"test2");
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
        let res = parse_sequences(
            seq(b"@test\nAGCT\n+test\n~~a!\n@test2\nTGCA\n+test\nWUI9"),
            |_| (),
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCT");
                        assert_eq!(&seq.qual.unwrap()[..], b"~~a!");
                    },
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"TGCA");
                        assert_eq!(&seq.qual.unwrap()[..], b"WUI9");
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(i, 2);
        assert_eq!(res, Ok(()));

        let mut i = 0;
        let res = parse_sequences(
            seq(b"@test\r\nAGCT\r\n+test\r\n~~a!\r\n@test2\r\nTGCA\r\n+test\r\nWUI9"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCT");
                        assert_eq!(&seq.qual.unwrap()[..], b"~~a!");
                    },
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"TGCA");
                        assert_eq!(&seq.qual.unwrap()[..], b"WUI9");
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
    fn test_fastq_endings() {
        //! Check for the absence of a panic. The parser previously assumed
        //! if the ID ended with an `\r\n` then the sequence did also.
        //! (Discovered via fuzzing)
        let res = parse_sequences(seq(b"@\r\n\n+A\n@"), |_| (), |_seq| {});
        assert_eq!(res, Ok(()));
    }

    #[test]
    fn test_wrapped_fasta() {
        let mut i = 0;
        let res = parse_sequences(
            seq(b">test\nAGCT\nTCG\n>test2\nG"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCTTCG");
                        assert_eq!(seq.qual, None);
                    },
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"G");
                        assert_eq!(seq.qual, None);
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(res, Ok(()));
        assert_eq!(i, 2);

        let mut i = 0;
        let res = parse_sequences(
            seq(b">test\r\nAGCT\r\nTCG\r\n>test2\r\nG"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCTTCG");
                        assert_eq!(seq.qual, None);
                    },
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"G");
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

    #[test]
    fn test_premature_endings() {
        let mut i = 0;
        let res = parse_sequences(
            seq(b">test\nAGCT\n>test2"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCT");
                        assert_eq!(seq.qual, None);
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(i, 1);
        let e = res.unwrap_err();
        assert_eq!(e.error_type, ParseErrorType::PrematureEOF);
        assert_eq!(e.record, 2);

        let mut i = 0;
        let res = parse_sequences(
            seq(b"@test\nAGCT\n+test\n~~a!\n@test2\nTGCA"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCT");
                        assert_eq!(&seq.qual.unwrap()[..], b"~~a!");
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(i, 1);
        let e = res.unwrap_err();
        assert_eq!(e.error_type, ParseErrorType::PrematureEOF);
        assert_eq!(e.record, 2);

        // we allow a few extra newlines at the ends of FASTQs
        let mut i = 0;
        let res = parse_sequences(
            seq(b"@test\nAGCT\n+test\n~~a!\n\n"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCT");
                        assert_eq!(&seq.qual.unwrap()[..], b"~~a!");
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(i, 1);
        assert_eq!(res, Ok(()));

        // but if there's additional data past the newlines it's an error
        let mut i = 0;
        let res = parse_sequences(
            seq(b"@test\nAGCT\n+test\n~~a!\n\n@TEST\nA\n+TEST\n~"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"AGCT");
                        assert_eq!(&seq.qual.unwrap()[..], b"~~a!");
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(i, 1);
        let e = res.unwrap_err();
        assert_eq!(e.error_type, ParseErrorType::PrematureEOF);
        assert_eq!(e.record, 2);

        // test that an abrupt stop in a FASTA triggers an error
        let mut i = 0;
        let res = parse_sequences(
            seq(b">test\nACGT\n>test2\n"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"test");
                        assert_eq!(&seq.seq[..], b"ACGT");
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(i, 1);
        let e = res.unwrap_err();
        assert_eq!(e.error_type, ParseErrorType::PrematureEOF);
        assert_eq!(e.record, 2);
    }

    #[test]
    fn test_empty_records() {
        let mut i = 0;
        let res = parse_sequences(
            seq(b"@\n\n+\n\n@test2\nTGCA\n+test2\n~~~~\n"),
            |stype| {
                assert_eq!(stype, "FASTQ");
            },
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"");
                        assert_eq!(&seq.seq[..], b"");
                        assert_eq!(&seq.qual.unwrap()[..], b"");
                    },
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"TGCA");
                        assert_eq!(&seq.qual.unwrap()[..], b"~~~~");
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(res, Ok(()));
        assert_eq!(i, 2);

        let mut i = 0;
        let res = parse_sequences(
            seq(b">\n\n>shine\nAGGAGGU"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"");
                        assert_eq!(&seq.seq[..], b"");
                        assert_eq!(seq.qual, None);
                    },
                    1 => {
                        assert_eq!(&seq.id[..], b"shine");
                        assert_eq!(&seq.seq[..], b"AGGAGGU");
                        assert_eq!(seq.qual, None);
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(i, 2);
        assert_eq!(res, Ok(()));

        let mut i = 0;
        let res = parse_sequences(
            seq(b">\r\n\r\n>shine\r\nAGGAGGU"),
            |_| {},
            |seq| {
                match i {
                    0 => {
                        assert_eq!(&seq.id[..], b"");
                        assert_eq!(&seq.seq[..], b"");
                        assert_eq!(seq.qual, None);
                    },
                    1 => {
                        assert_eq!(&seq.id[..], b"shine");
                        assert_eq!(&seq.seq[..], b"AGGAGGU");
                        assert_eq!(seq.qual, None);
                    },
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(i, 2);
        assert_eq!(res, Ok(()));
    }

    #[test]
    fn test_buffer() {
        let mut buf: RecBuffer<FASTA> = RecBuffer::from_bytes(b">test\nACGT");
        let rec = buf.next().unwrap().unwrap();
        assert_eq!(rec.id, b"test", "Record has the right ID");
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
            let mut rec_buffer = rec_reader.get_buffer::<FASTQ>(0);
            for _s in rec_buffer.by_ref() {
                // record is incomplete
                panic!("No initial record should be parsed")
            }
            rec_buffer.pos
        };

        // refill the buffer, but we're not done quite yet
        assert_eq!(rec_reader.refill(used).unwrap(), false);

        // now we should see both records
        let mut rec_buffer = rec_reader.get_buffer::<FASTQ>(0);

        // there should be a record assuming the parser
        // handled the buffer boundary
        let iterated_seq = rec_buffer.by_ref().next();
        let seq = iterated_seq.unwrap();
        assert_eq!(seq.unwrap().id, b"A");

        // but not another because the buffer's too short
        let iterated_seq = rec_buffer.by_ref().next();
        assert!(iterated_seq.is_none());

        // TODO: refill and check for the last record
    }
}
