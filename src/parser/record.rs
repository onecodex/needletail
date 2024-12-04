use std::borrow::Cow;
use std::io::Write;

use memchr::memchr;

use crate::errors::ParseError;
use crate::parser::fasta::BufferPosition as FastaBufferPosition;
use crate::parser::fastq::BufferPosition as FastqBufferPosition;
use crate::parser::utils::{Format, LineEnding, Position};
use crate::Sequence;

#[derive(Debug, Clone)]
enum BufferPositionKind<'a> {
    Fasta(&'a FastaBufferPosition),
    Fastq(&'a FastqBufferPosition),
}

/// A FASTA or FASTQ record
#[derive(Debug, Clone)]
pub struct SequenceRecord<'a> {
    buffer: &'a [u8],
    buf_pos: BufferPositionKind<'a>,
    position: &'a Position,
    line_ending: LineEnding,
}

impl<'a> SequenceRecord<'a> {
    pub(crate) fn new_fasta(
        buffer: &'a [u8],
        buf_pos: &'a FastaBufferPosition,
        position: &'a Position,
        line_ending: Option<LineEnding>,
    ) -> Self {
        Self {
            buffer,
            position,
            buf_pos: BufferPositionKind::Fasta(buf_pos),
            line_ending: line_ending.unwrap_or(LineEnding::Unix),
        }
    }

    pub(crate) fn new_fastq(
        buffer: &'a [u8],
        buf_pos: &'a FastqBufferPosition,
        position: &'a Position,
        line_ending: Option<LineEnding>,
    ) -> Self {
        Self {
            buffer,
            position,
            buf_pos: BufferPositionKind::Fastq(buf_pos),
            line_ending: line_ending.unwrap_or(LineEnding::Unix),
        }
    }

    /// Returns the format of the record
    #[inline]
    pub fn format(&self) -> Format {
        match self.buf_pos {
            BufferPositionKind::Fasta(_) => Format::Fasta,
            BufferPositionKind::Fastq(_) => Format::Fastq,
        }
    }

    /// Returns the id of the record
    #[inline]
    pub fn id(&self) -> &[u8] {
        match self.buf_pos {
            BufferPositionKind::Fasta(bp) => bp.id(self.buffer),
            BufferPositionKind::Fastq(bp) => bp.id(self.buffer),
        }
    }

    /// Returns the raw sequence of the record. Only matters for FASTA since it can contain
    /// newlines.
    #[inline]
    pub fn raw_seq(&self) -> &[u8] {
        match self.buf_pos {
            BufferPositionKind::Fasta(bp) => bp.raw_seq(self.buffer),
            BufferPositionKind::Fastq(bp) => bp.seq(self.buffer),
        }
    }

    /// Returns the cleaned up sequence of the record. For FASTQ it is the same as `raw_seq` but
    /// for FASTA it is `raw_seq` minus all the `\r\n`
    pub fn seq(&self) -> Cow<[u8]> {
        match self.buf_pos {
            BufferPositionKind::Fasta(bp) => bp.seq(self.buffer),
            BufferPositionKind::Fastq(bp) => bp.seq(self.buffer).into(),
        }
    }

    /// Returns the quality line if there is one.
    /// Always `None` for FASTA and `Some` for FASTQ, even if the quality line is empty.
    #[inline]
    pub fn qual(&self) -> Option<&[u8]> {
        match self.buf_pos {
            BufferPositionKind::Fasta(_) => None,
            BufferPositionKind::Fastq(bp) => Some(bp.qual(self.buffer)),
        }
    }

    /// Returns the full sequence, including line endings. This doesn't include a trailing newline.
    #[inline]
    pub fn all(&self) -> &[u8] {
        match self.buf_pos {
            BufferPositionKind::Fasta(bp) => bp.all(self.buffer),
            BufferPositionKind::Fastq(bp) => bp.all(self.buffer),
        }
    }

    /// Return the number of bases in the sequence, computed efficiently.
    #[inline]
    pub fn num_bases(&self) -> usize {
        match self.buf_pos {
            BufferPositionKind::Fasta(bp) => bp.num_bases(self.buffer),
            BufferPositionKind::Fastq(bp) => bp.num_bases(self.buffer),
        }
    }

    /// Return the line number in the file of the start of the sequence
    pub fn start_line_number(&self) -> u64 {
        self.position.line
    }

    /// Return the line/byte position of the start of the sequence
    pub fn position(&self) -> &Position {
        self.position
    }

    /// Which line ending is this record using?
    pub fn line_ending(&self) -> LineEnding {
        self.line_ending
    }

    /// Write record back to a `Write` instance. By default it will use the original line ending but
    /// you can force it to use another one.
    pub fn write(
        &self,
        writer: &mut dyn Write,
        forced_line_ending: Option<LineEnding>,
    ) -> Result<(), ParseError> {
        match self.buf_pos {
            BufferPositionKind::Fasta(_) => write_fasta(
                self.id(),
                self.raw_seq(),
                writer,
                forced_line_ending.unwrap_or(self.line_ending),
            ),
            BufferPositionKind::Fastq(_) => write_fastq(
                self.id(),
                self.raw_seq(),
                self.qual(),
                writer,
                forced_line_ending.unwrap_or(self.line_ending),
            ),
        }
    }
}

impl<'a> Sequence<'a> for SequenceRecord<'a> {
    fn sequence(&'a self) -> &'a [u8] {
        self.raw_seq()
    }
}

/// Mask tabs in header lines to `|`s
pub fn mask_header_tabs(id: &[u8]) -> Option<Vec<u8>> {
    memchr(b'\t', id).map(|_| {
        id.iter()
            .map(|x| if *x == b'\t' { b'|' } else { *x })
            .collect()
    })
}

/// Convert bad UTF8 characters into �s
pub fn mask_header_utf8(id: &[u8]) -> Option<Vec<u8>> {
    // this may potentially change the length of the id; we should probably
    // be doing something trickier like converting
    match String::from_utf8_lossy(id) {
        Cow::Owned(s) => Some(s.into_bytes()),
        Cow::Borrowed(_) => None,
    }
}

/// Write a FASTA record
pub fn write_fasta(
    id: &[u8],
    seq: &[u8],
    writer: &mut dyn Write,
    line_ending: LineEnding,
) -> Result<(), ParseError> {
    let ending = line_ending.to_bytes();
    writer.write_all(b">")?;
    writer.write_all(id)?;
    writer.write_all(&ending)?;
    writer.write_all(seq)?;
    writer.write_all(&ending)?;
    Ok(())
}

pub fn write_fastq(
    id: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
    writer: &mut dyn Write,
    line_ending: LineEnding,
) -> Result<(), ParseError> {
    let ending = line_ending.to_bytes();
    writer.write_all(b"@")?;
    writer.write_all(id)?;
    writer.write_all(&ending)?;
    writer.write_all(seq)?;
    writer.write_all(&ending)?;
    writer.write_all(b"+")?;
    writer.write_all(&ending)?;
    // this is kind of a hack, but we want to allow writing out sequences
    // that don't have qualitys so this will mask to "good" if the quality
    // slice is empty
    if let Some(qual) = qual {
        writer.write_all(qual)?;
    } else {
        writer.write_all(&vec![b'I'; seq.len()])?;
    }
    writer.write_all(&ending)?;
    Ok(())
}

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use crate::parse_fastx_reader;
    fn seq(s: &[u8]) -> Cursor<&[u8]> {
        Cursor::new(s)
    }

    #[test]
    fn test_start_line_number() {
        let mut reader =
            parse_fastx_reader(seq(b"@test\nACGT\n+\nIIII\n@test2\nACGT\n+\nIIII")).unwrap();

        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.start_line_number(), 1);

        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.start_line_number(), 5);
    }

    #[test]
    fn test_position() {
        let mut reader = parse_fastx_reader(seq(
            b"@test1\nACGT\n+\nIIII\n@test222\nACGT\n+\nIIII\n@test3\nACGT\n+\nIIII",
        ))
        .unwrap();

        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.position().byte(), 0);

        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.position().byte(), 19);

        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.position().byte(), 40);
    }
}
