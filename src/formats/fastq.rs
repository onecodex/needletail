use std::borrow::Cow;
use std::cmp::min;
use std::io::Write;

use memchr::memchr;

use crate::buffer::RecBuffer;
use crate::seq::Sequence;
use crate::util::{memchr_both, ParseError, ParseErrorType};


#[derive(Debug)]
pub struct FASTQ<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub id2: &'a [u8],
    pub qual: &'a [u8],
}

impl<'a> FASTQ<'a> {
    pub fn write<W>(&self, mut writer: W) -> Result<(), ParseError> where W: Write {
        writer.write(b"@")?;
        writer.write(&self.id)?;
        writer.write(b"\n")?;
        writer.write(&self.seq)?;
        writer.write(b"+\n")?;
        if self.seq.len() != self.qual.len() {
            writer.write(&vec![b'I'; self.seq.len()])?;
        } else {
            writer.write(&self.qual)?;
        }
        writer.write(b"\n")?;
        Ok(())
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
                let context = String::from_utf8_lossy(&buf[0..min(16, buf.len())]);
                let e =
                    ParseError::new("Record must start with '@'", ParseErrorType::InvalidHeader)
                        .record(self.count)
                        .context(context);
                return Some(Err(e));
            }
        }

        let id_end;
        match memchr(b'\n', &buf) {
            Some(i) => id_end = i + 1,
            None => return None,
        };
        let mut id = &buf[1..id_end - 1];

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
        let id2 = &buf[seq_end..id2_end - 1];

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
            if qual_end != buf.len() + 1 + windows_ending {
                return None;
            }
            buffer_used -= 1 + windows_ending;
            qual_end -= windows_ending;
        }
        let mut qual = &buf[id2_end..qual_end - 1];

        // clean up any extra '\r' from the id and seq
        if !id.is_empty() && id[id.len() - 1] == b'\r' {
            id = &id[..id.len() - 1];
        }
        if !seq.is_empty() && seq[seq.len() - 1] == b'\r' {
            seq = &seq[..seq.len() - 1];
        }
        // we do qual separately in case this is the end of the file
        if !qual.is_empty() && qual[qual.len() - 1] == b'\r' {
            qual = &qual[..qual.len() - 1];
        }

        self.pos += buffer_used;
        self.count += 1;
        Some(Ok(FASTQ { id, seq, id2, qual }))
    }
}

impl<'a> From<FASTQ<'a>> for Sequence<'a> {
    fn from(fastq: FASTQ<'a>) -> Sequence<'a> {
        let qual = if fastq.seq.len() != fastq.qual.len() {
            None
        } else {
            Some(fastq.qual)
        };
        Sequence::new(fastq.id, Cow::from(fastq.seq), qual)
    }
}

impl<'a> From<&'a Sequence<'a>> for FASTQ<'a> {
    fn from(seq: &'a Sequence<'a>) -> FASTQ<'a> {
        let qual = match &seq.qual {
            None => &b""[..],
            Some(q) => &q,
        };
        FASTQ {
            id: &seq.id,
            seq: &seq.seq,
            id2: b"",
            qual: qual,
        }
    }
}
