use std::cmp::min;
use std::io::Write;

use memchr::memchr;

use crate::formats::buffer::RecParser;
use crate::seq::{Sequence, SequenceRecord};
use crate::util::{memchr_both_last, ParseError, ParseErrorType};

#[derive(Debug)]
pub struct FastaRecord<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
}

impl<'a> FastaRecord<'a> {
    pub fn write(&self, writer: &mut dyn Write, ending: &[u8]) -> Result<(), ParseError> {
        writer.write_all(b">")?;
        writer.write_all(&self.id)?;
        writer.write_all(ending)?;
        writer.write_all(&self.seq)?;
        writer.write_all(ending)?;
        Ok(())
    }
}

impl<'a> Sequence<'a> for FastaRecord<'a> {
    fn sequence(&self) -> &'a [u8] {
        self.seq
    }
}

impl<'a> From<FastaRecord<'a>> for SequenceRecord<'a> {
    fn from(fasta: FastaRecord<'a>) -> SequenceRecord<'a> {
        SequenceRecord::new(fasta.id.into(), fasta.seq.into(), None)
    }
}

pub struct FastaParser<'a> {
    buf: &'a [u8],
    last: bool,
    pos: usize,
}

impl<'a> FastaParser<'a> {
    pub fn new(buf: &'a [u8], last: bool) -> Result<Self, ParseError> {
        if buf[0] != b'>' {
            let context = String::from_utf8_lossy(&buf[..min(32, buf.len())]);
            return Err(ParseError::new(
                "FASTA record must start with '>'",
                ParseErrorType::InvalidHeader,
            )
            .context(context));
        }

        Ok(FastaParser { buf, last, pos: 0 })
    }
}

impl<'a> Iterator for FastaParser<'a> {
    type Item = Result<FastaRecord<'a>, ParseError>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let buf = &self.buf[self.pos..];
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
        match (memchr_both_last(b'\n', b'>', &buf[id_end..]), self.last) {
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
            .context(context)));
        }
        let mut seq = &buf[id_end..seq_end];
        if seq.len() > 0 && seq[seq.len() - 1] == b'\n' {
            seq = &seq[..seq.len() - 1];
        }
        if seq.len() > 0 && seq[seq.len() - 1] == b'\r' {
            seq = &seq[..seq.len() - 1];
        }

        self.pos += seq_end;
        Some(Ok(FastaRecord { id, seq }))
    }
}

impl<'a> RecParser<'a> for FastaParser<'a> {
    type Header = ();

    fn from_buffer(buf: &[u8], last: bool) -> FastaParser {
        FastaParser { buf, last, pos: 0 }
    }

    fn header(&mut self) -> Result<Self::Header, ParseError> {
        Ok(())
    }

    fn eof(&self) -> Result<(), ParseError> {
        check_end(&self.buf[self.pos..], self.last)
    }

    fn used(&self) -> usize {
        self.pos
    }
}

pub fn check_end(buf: &[u8], last: bool) -> Result<(), ParseError> {
    // check if there's anything left stuff in the buffer (besides returns)
    if !last {
        return Err(
            ParseError::new("File ended abruptly", ParseErrorType::PrematureEOF),
            // .record(count + 1),
        );
    }
    for c in &buf[..] {
        if c != &b'\r' && c != &b'\n' {
            let end = min(32, buf.len());
            let context = String::from_utf8_lossy(&buf[..end]);
            return Err(ParseError::new(
                "File had extra data past end of records",
                ParseErrorType::PrematureEOF,
            )
            // .record(count + 1)
            .context(context));
        }
    }
    Ok(())
}

#[cfg(test)]
mod test {
    use std::fs::File;
    use std::io::Cursor;
    use std::path::Path;

    use super::FastaParser;
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
                    }
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"GATC");
                        assert_eq!(seq.qual, None);
                    }
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
                    }
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"TAGC");
                        assert_eq!(seq.qual, None);
                    }
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
        assert_eq!(e.msg, String::from("Could not detect file type"));
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
                        }
                        1 => {
                            assert_eq!(&seq.id[..], b"test2");
                            assert_eq!(&seq.seq[..], b"TAGC");
                            assert_eq!(seq.qual, None);
                        }
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
                    }
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"G");
                        assert_eq!(seq.qual, None);
                    }
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
                    }
                    1 => {
                        assert_eq!(&seq.id[..], b"test2");
                        assert_eq!(&seq.seq[..], b"G");
                        assert_eq!(seq.qual, None);
                    }
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
                    }
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
                    }
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
                    }
                    1 => {
                        assert_eq!(&seq.id[..], b"shine");
                        assert_eq!(&seq.seq[..], b"AGGAGGU");
                        assert_eq!(seq.qual, None);
                    }
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
                    }
                    1 => {
                        assert_eq!(&seq.id[..], b"shine");
                        assert_eq!(&seq.seq[..], b"AGGAGGU");
                        assert_eq!(seq.qual, None);
                    }
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            },
        );
        assert_eq!(i, 2);
        assert_eq!(res, Ok(()));
    }

    #[test]
    fn test_reader() {
        let mut reader = FastaParser::new(b">test\nACGT", true).unwrap();
        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.id, b"test", "Record has the right ID");
        assert_eq!(rec.seq, b"ACGT", "Record has the right sequence");

        let mut reader = FastaParser::new(b">test", true).unwrap();
        assert!(reader.next().is_none(), "Incomplete record returns None");
    }
}
