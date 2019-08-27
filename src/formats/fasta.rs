use std::io::Write;

use memchr::memchr;

use crate::buffer::RecBuffer;
use crate::seq::Sequence;
use crate::util::{memchr_both, strip_whitespace, ParseError, ParseErrorType};

#[derive(Debug)]
pub struct FASTA<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
}

impl<'a> FASTA<'a> {
    pub fn write<W>(&self, mut writer: W) -> Result<(), ParseError> where W: Write {
        writer.write(b">")?;
        writer.write(&self.id)?;
        writer.write(b"\n")?;
        writer.write(&self.seq)?;
        writer.write(b"\n")?;
        Ok(())
    }
}

impl<'a> Iterator for RecBuffer<'a, FASTA<'static>> {
    type Item = Result<FASTA<'a>, ParseError>;

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
        match (memchr_both(b'\n', b'>', &buf[id_end..]), self.last) {
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
            .record(self.count + 1)
            .context(context)));
        }
        let mut seq = &buf[id_end..seq_end];
        if seq[seq.len() - 1] == b'\r' {
            seq = &seq[..seq.len()];
        }

        self.pos += seq_end;
        self.count += 1;
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
