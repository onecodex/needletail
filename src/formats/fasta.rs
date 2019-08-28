use std::io::Write;

use memchr::memchr;

use crate::buffer::{RecBuffer, RecordFormat};
use crate::seq::Sequence;
use crate::util::{memchr_both, strip_whitespace, ParseError, ParseErrorType};

#[derive(Debug)]
pub struct FASTA<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
}

impl<'a> FASTA<'a> {
    pub fn write<W>(&self, mut writer: W) -> Result<(), ParseError>
    where
        W: Write,
    {
        writer.write(b">")?;
        writer.write(&self.id)?;
        writer.write(b"\n")?;
        writer.write(&self.seq)?;
        writer.write(b"\n")?;
        Ok(())
    }
}

impl<'a> RecordFormat<'a> for FASTA<'a> {
    fn parse(rbuf: &'a mut RecBuffer) -> Option<Result<Self, ParseError>> {
        let buf = &rbuf.buf[rbuf.pos..];
        if buf.is_empty() {
            return None;
        }

        let id_end;
        match memchr(b'\n', &buf) {
            Some(i) => id_end = i + 1,
            None => return None,
        };
        let mut id = &buf[1..id_end - 1];
        if !id.is_empty() && id[id.len() - 1] == b'\r' {
            id = &id[..id.len() - 1];
        }

        let seq_end;
        match (memchr_both(b'\n', b'>', &buf[id_end..]), rbuf.last) {
            (Some(i), _) => seq_end = id_end + i + 1,
            (None, true) => seq_end = buf.len(),
            (None, false) => return None,
        };
        if id_end == seq_end {
            let context = String::from_utf8_lossy(id);
            return Some(Err(ParseError::new(
                "Sequence completely empty",
                ParseErrorType::PrematureEOF,
            )
            .record(rbuf.count + 1)
            .context(context)));
        }
        let mut seq = &buf[id_end..seq_end];
        if seq[seq.len() - 1] == b'\r' {
            seq = &seq[..seq.len()];
        }

        rbuf.pos += seq_end;
        rbuf.count += 1;
        Some(Ok(FASTA { id, seq }))
    }
}

impl<'a> From<FASTA<'a>> for Sequence<'a> {
    fn from(fasta: FASTA<'a>) -> Sequence<'a> {
        Sequence::new(fasta.id, strip_whitespace(fasta.seq), None)
    }
}

impl<'a> From<&'a Sequence<'a>> for FASTA<'a> {
    fn from(seq: &'a Sequence<'a>) -> FASTA<'a> {
        FASTA {
            id: &seq.id,
            seq: &seq.seq,
        }
    }
}

#[cfg(test)]
mod test {
    use std::fs::File;
    use std::io::Cursor;
    use std::path::Path;

    use super::FASTA;
    use crate::buffer::RecBuffer;
    use crate::formats::parse_sequences;
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
}
