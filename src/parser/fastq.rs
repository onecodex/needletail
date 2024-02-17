//! The vast majority of the code is taken from https://github.com/markschl/seq_io/blob/master/src/fastq.rs

use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

use crate::errors::{ErrorPosition, ParseError};
use crate::parser::record::SequenceRecord;
use crate::parser::utils::{
    fill_buf, find_line_ending, grow_to, trim_cr, FastxReader, Format, LineEnding, Position,
    BUFSIZE,
};
use memchr::memchr;

/// Represents the position of a record within a buffer
#[derive(Debug, Clone, Default)]
pub struct BufferPosition {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) seq: usize,
    pub(crate) sep: usize,
    pub(crate) qual: usize,
}

impl BufferPosition {
    #[inline]
    pub(crate) fn is_new(&self) -> bool {
        self.end == 0
    }

    #[inline]
    pub(crate) fn len(&self) -> u64 {
        (self.end + 1 - self.start) as u64
    }

    #[inline]
    pub(crate) fn id<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.start + 1..self.seq - 1])
    }

    #[inline]
    pub(crate) fn seq<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.seq..self.sep - 1])
    }

    #[inline]
    pub(crate) fn qual<'a>(&'a self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.qual..self.end])
    }

    #[inline]
    pub(crate) fn num_bases<'a>(&'a self, buffer: &'a [u8]) -> usize {
        self.seq(buffer).len()
    }

    #[inline]
    fn find_line_ending<'a>(&'a self, buffer: &'a [u8]) -> Option<LineEnding> {
        find_line_ending(self.all(buffer))
    }

    #[inline]
    pub(crate) fn all<'a>(&self, buffer: &'a [u8]) -> &'a [u8] {
        &buffer[self.start..self.end]
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Ord, PartialOrd)]
enum SearchPosition {
    Id,
    Sequence,
    Separator,
    Quality,
}

/// Parser for FASTQ files.
/// Only use this directly if you know your file is FASTQ and that it is not compressed as
/// it does not handle decompression.
/// If you are unsure, it's better to use [`parse_fastx_file`](fn.parse_fastx_file.html).
pub struct Reader<R: io::Read> {
    buf_reader: buffer_redux::BufReader<R>,
    buf_pos: BufferPosition,
    search_pos: SearchPosition,
    position: Position,
    finished: bool,
    line_ending: Option<LineEnding>,
}

impl<R> Reader<R>
where
    R: io::Read,
{
    /// Creates a new reader with the default buffer size of 64 KiB
    ///
    /// # Example:
    ///
    /// ```
    /// use needletail::parser::{FastqReader, FastxReader};
    /// let fastq = b"@id\nACGT\n+\nIIII";
    ///
    /// let mut reader = FastqReader::new(&fastq[..]);
    /// let record = reader.next().unwrap().unwrap();
    /// assert_eq!(record.id(), b"id")
    /// ```
    pub fn new(reader: R) -> Self {
        Self::with_capacity(reader, BUFSIZE)
    }

    /// Creates a new reader with a given buffer capacity. The minimum allowed
    /// capacity is 3.
    pub fn with_capacity(reader: R, capacity: usize) -> Self {
        assert!(capacity >= 3);
        Self {
            buf_reader: buffer_redux::BufReader::with_capacity(capacity, reader),
            buf_pos: BufferPosition::default(),
            search_pos: SearchPosition::Id,
            position: Position::new(1, 0),
            finished: false,
            line_ending: None,
        }
    }
}

impl Reader<File> {
    /// Creates a reader from a file path.
    ///
    /// # Example:
    ///
    /// ```no_run
    /// use needletail::parser::{FastxReader, FastqReader};
    ///
    /// let mut reader = FastqReader::from_path("seqs.fastq").unwrap();
    ///
    /// // (... do something with the reader)
    /// ```
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        File::open(path).map(Self::new)
    }
}

impl<R> Reader<R>
where
    R: io::Read,
{
    #[inline]
    fn get_buf(&self) -> &[u8] {
        self.buf_reader.buffer()
    }

    // TODO: avoid duplication with find_incomplete.
    // TODO: having a single fn and adding branches introduce a noticeable slowdown
    /// Reads the current record and returns true if found.
    /// Returns false if incomplete because end of buffer reached,
    /// meaning that the last record may be incomplete.
    /// Updates `self.search_pos`.
    fn find(&mut self) -> Result<bool, ParseError> {
        self.buf_pos.seq = if let Some(p) = self.find_line(self.buf_pos.start) {
            p
        } else {
            self.search_pos = SearchPosition::Id;
            return Ok(false);
        };

        self.buf_pos.sep = if let Some(p) = self.find_line(self.buf_pos.seq) {
            p
        } else {
            self.search_pos = SearchPosition::Sequence;
            return Ok(false);
        };

        self.buf_pos.qual = if let Some(p) = self.find_line(self.buf_pos.sep) {
            p
        } else {
            self.search_pos = SearchPosition::Separator;
            return Ok(false);
        };

        self.buf_pos.end = if let Some(p) = self.find_line(self.buf_pos.qual) {
            p - 1
        } else {
            self.search_pos = SearchPosition::Quality;
            return Ok(false);
        };

        self.validate()?;

        Ok(true)
    }

    // Resumes reading an incomplete record without
    // re-searching positions that were already found.
    // The resulting position may still be incomplete (-> false).
    fn find_incomplete(&mut self) -> Result<bool, ParseError> {
        if self.search_pos == SearchPosition::Id {
            self.buf_pos.seq = if let Some(p) = self.find_line(self.buf_pos.start) {
                p
            } else {
                self.search_pos = SearchPosition::Id;
                return Ok(false);
            };
        }

        if self.search_pos <= SearchPosition::Sequence {
            self.buf_pos.sep = if let Some(p) = self.find_line(self.buf_pos.seq) {
                p
            } else {
                self.search_pos = SearchPosition::Sequence;
                return Ok(false);
            };
        }

        if self.search_pos <= SearchPosition::Separator {
            self.buf_pos.qual = if let Some(p) = self.find_line(self.buf_pos.sep) {
                p
            } else {
                self.search_pos = SearchPosition::Separator;
                return Ok(false);
            };
        }

        if self.search_pos <= SearchPosition::Quality {
            self.buf_pos.end = if let Some(p) = self.find_line(self.buf_pos.qual) {
                p - 1
            } else {
                self.search_pos = SearchPosition::Quality;
                return Ok(false);
            };
        }

        self.search_pos = SearchPosition::Id;

        self.validate()?;

        Ok(true)
    }

    /// Verify that the record is valid:
    /// - starts with @
    /// - separator line starts with -
    /// - quality and sequence have the same length
    fn validate(&mut self) -> Result<(), ParseError> {
        let start_byte = self.get_buf()[self.buf_pos.start];
        if start_byte != b'@' {
            self.finished = true;
            return Err(ParseError::new_invalid_start(
                start_byte,
                self.get_error_pos(0, false),
                Format::Fastq,
            ));
        }

        let sep_byte = self.get_buf()[self.buf_pos.sep];
        if sep_byte != b'+' {
            self.finished = true;
            return Err(ParseError::new_invalid_separator(
                sep_byte,
                self.get_error_pos(2, true),
            ));
        }

        let buf = self.get_buf();
        // We assume we only have ASCII in sequence and quality
        let seq_len = self.buf_pos.seq(buf).len();
        let qual_len = self.buf_pos.qual(buf).len();

        // TODO: we don't do that every time because it's a ~90% performance penalty.
        // TODO: mention it on the README
        // And we can further validate quality chars
        // and the vast majority of files don't have this issue
        // let qual_len = self
        //     .buf_pos
        //     .qual(&buf)
        //     .iter()
        //     .filter(|c| *c >= &b'!' && *c <= &b'~')
        //     .count();

        if seq_len != qual_len {
            self.finished = true;
            return Err(ParseError::new_unequal_length(
                seq_len,
                qual_len,
                self.get_error_pos(0, true),
            ));
        }
        Ok(())
    }

    fn get_error_pos(&self, line_offset: u64, parse_id: bool) -> ErrorPosition {
        let id = if parse_id && self.buf_pos.seq - self.buf_pos.start > 1 {
            let id = self
                .buf_pos
                .id(self.get_buf())
                .split(|b| *b == b' ')
                .next()
                .unwrap();
            Some(String::from_utf8_lossy(id).into())
        } else {
            None
        };
        ErrorPosition {
            line: self.position.line() + line_offset,
            id,
        }
    }

    #[inline]
    fn find_line(&self, search_start: usize) -> Option<usize> {
        memchr(b'\n', &self.get_buf()[search_start..]).map(|pos| search_start + pos + 1)
    }

    /// Called when we couldn't find a complete record.
    /// We might be at EOF, buffer might be too small or we need to refill it
    fn next_complete(&mut self) -> Result<bool, ParseError> {
        loop {
            if self.get_buf().len() < self.buf_reader.capacity() {
                // EOF reached, there will be no next record
                return self.check_end();
            }

            if self.buf_pos.start == 0 {
                // first record already incomplete -> buffer too small
                self.grow();
            } else {
                // not the first record -> buffer may be big enough but we need to make some space
                self.make_room();
            }

            fill_buf(&mut self.buf_reader)?;

            if self.find_incomplete()? {
                return Ok(true);
            }
        }
    }

    /// Checks for EOF.
    /// If there is one last record that can be sent, return `true` otherwise `false`.
    fn check_end(&mut self) -> Result<bool, ParseError> {
        self.finished = true;
        if self.search_pos == SearchPosition::Quality {
            // no line ending at end of last record
            self.buf_pos.end = self.get_buf().len();
            self.validate()?;
            return Ok(true);
        }

        // It allows some blank lines at the end of the file
        let rest = &self.get_buf()[self.buf_pos.start..];
        if rest.split(|c| *c == b'\n').all(|l| trim_cr(l).is_empty()) {
            return Ok(false);
        }

        Err(ParseError::new_unexpected_end(
            self.get_error_pos(self.search_pos as u64, self.search_pos > SearchPosition::Id),
            Format::Fastq,
        ))
    }

    // Grow the internal buffer. Used if the original buffer is not big
    // enough for a record
    fn grow(&mut self) {
        let cap = self.buf_reader.capacity();
        let new_size = grow_to(cap);
        let additional = new_size - cap;
        self.buf_reader.reserve(additional);
    }

    // Consume bytes from records we've seen and move incomplete bytes to start of buffer
    fn make_room(&mut self) {
        let consumed = self.buf_pos.start;
        self.buf_reader.consume(consumed);
        self.buf_reader.make_room();

        self.buf_pos.start = 0;

        if self.search_pos >= SearchPosition::Sequence {
            self.buf_pos.seq -= consumed;
        }
        if self.search_pos >= SearchPosition::Separator {
            self.buf_pos.sep -= consumed;
        }
        if self.search_pos >= SearchPosition::Quality {
            self.buf_pos.qual -= consumed;
        }
    }
}

impl<R: io::Read + Send> FastxReader for Reader<R> {
    fn next(&mut self) -> Option<Result<SequenceRecord, ParseError>> {
        // No more records to read
        if self.finished {
            return None;
        }

        // Empty buffer, let's fill it
        if self.get_buf().is_empty() {
            // If we get an ParseError when reading or get back 0 bytes, we're done
            match fill_buf(&mut self.buf_reader) {
                Ok(n) => {
                    if n == 0 {
                        self.finished = true;
                        return None;
                    }
                }
                Err(e) => {
                    return Some(Err(e.into()));
                }
            };
        }

        // If we already did look at a record, let's setup for the next one
        if !self.buf_pos.is_new() {
            self.position.byte += self.buf_pos.len();
            self.position.line += 4;
            self.buf_pos.start = self.buf_pos.end + 1;
        }

        // Can we identify all the positions of each element of the next record?
        let complete = match self.find() {
            Ok(f) => f,
            Err(e) => {
                return Some(Err(e));
            }
        };

        // If it's not complete, try to fetch more from the buffer until we have it in full
        if !complete {
            // Did we get a record?
            let got_record = match self.next_complete() {
                Ok(f) => f,
                Err(e) => {
                    return Some(Err(e));
                }
            };

            if !got_record {
                return None;
            }
        }
        if self.line_ending.is_none() {
            self.line_ending = self.buf_pos.find_line_ending(self.get_buf());
        }
        // We got one!
        Some(Ok(SequenceRecord::new_fastq(
            self.get_buf(),
            &self.buf_pos,
            &self.position,
            self.line_ending,
        )))
    }

    fn position(&self) -> &Position {
        &self.position
    }

    fn line_ending(&self) -> Option<LineEnding> {
        self.line_ending
    }
}

#[cfg(test)]
mod test {
    use std::io::Cursor;

    use super::Reader;
    use crate::errors::ParseErrorKind;
    use crate::parser::utils::LineEnding;
    use crate::FastxReader;

    fn seq(s: &[u8]) -> Cursor<&[u8]> {
        Cursor::new(s)
    }

    #[test]
    fn test_simple_fastq() {
        // Test both line endings
        let sequences = vec![
            (
                "@test\nAGCT\n+test\n~~a!\n@test2\nTGCA\n+test\nWUI9",
                LineEnding::Unix,
            ),
            (
                "@test\r\nAGCT\r\n+test\r\n~~a!\r\n@test2\r\nTGCA\r\n+test\r\nWUI9",
                LineEnding::Windows,
            ),
        ];

        for (sequence, line_ending) in sequences {
            let mut i = 0;
            let mut reader = Reader::new(seq(sequence.as_bytes()));
            while let Some(record) = reader.next() {
                let rec = record.unwrap();
                match i {
                    0 => {
                        assert_eq!(&rec.id(), b"test");
                        assert_eq!(&rec.raw_seq(), b"AGCT");
                        assert_eq!(&rec.qual().unwrap(), b"~~a!");
                        assert_eq!(reader.line_ending().unwrap(), line_ending);
                    }
                    1 => {
                        assert_eq!(&rec.id(), b"test2");
                        assert_eq!(&rec.raw_seq(), b"TGCA");
                        assert_eq!(&rec.qual().unwrap(), b"WUI9");
                        assert_eq!(reader.line_ending().unwrap(), line_ending);
                    }
                    _ => unreachable!("Too many records"),
                }
                i += 1;
            }
            assert_eq!(i, 2);
        }
    }

    #[test]
    fn test_eof_in_qual() {
        let mut reader = Reader::new(seq(b"@test\nACGT\n+\nIII"));
        let rec = reader.next().unwrap();
        assert!(rec.is_err());
        let e = rec.unwrap_err();
        // Not a eof error due to the way the validate is implemented
        assert_eq!(e.kind, ParseErrorKind::UnequalLengths);
    }

    #[test]
    fn test_eof_in_seq() {
        let mut reader = Reader::new(seq(b"@test\nAGCT\n+test\n~~a!\n@test2\nTGCA"));
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let rec2 = reader.next().unwrap();
        assert!(rec2.is_err());
        let e = rec2.unwrap_err();
        assert_eq!(e.kind, ParseErrorKind::UnexpectedEnd);
    }

    #[test]
    fn test_extra_empty_newlines_at_end_are_ok() {
        let mut reader = Reader::new(seq(b"@test\nAGCT\n+test\n~~a!\n\n"));
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_extra_non_empty_newlines_at_end_are_not_ok() {
        let mut reader = Reader::new(seq(b"@test\nAGCT\n+test\n~~a!\n\n@TEST\nA\n+TEST\n~"));
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let rec2 = reader.next().unwrap();
        let e = rec2.unwrap_err();
        assert_eq!(e.kind, ParseErrorKind::InvalidStart);
    }

    #[test]
    fn test_empty_records() {
        let mut reader = Reader::new(seq(b"@\n\n+\n\n@test2\nTGCA\n+test2\n~~~~\n"));
        let mut i = 0;
        while let Some(record) = reader.next() {
            let rec = record.unwrap();
            match i {
                0 => {
                    assert_eq!(&rec.id(), b"");
                    assert_eq!(&rec.raw_seq(), b"");
                    assert_eq!(&rec.qual().unwrap(), b"");
                    assert_eq!(rec.all(), b"@\n\n+\n");
                }
                1 => {
                    assert_eq!(&rec.id(), b"test2");
                    assert_eq!(&rec.raw_seq(), b"TGCA");
                    assert_eq!(&rec.qual().unwrap(), b"~~~~");
                    assert_eq!(rec.all(), b"@test2\nTGCA\n+test2\n~~~~");
                }
                _ => unreachable!("Too many records"),
            }
            i += 1
        }
        assert_eq!(i, 2);
    }

    #[test]
    fn test_weird_ncbi_file() {
        let test = b"@NCBI actually has files like this\nACGTACGATCGTACGTAGCTGCTAGCTAGCATGCATGACACACACGTACGATCGTACGTAGCTGCTAGCTAGCATGCATGACACAC\n+\n00000000000000000000000000000000000000000000000000000000000000000000000000000000000000\n@NCBI actually has files like this\n\n+\n\n@NCBI actually has files like this\nACGTACGATCGTACGTAGCTGCTAGCTAGCATGCATGACACACACGTACGATCGTACGTAGCTGCTAGCTAGCATGCATGACACAC\n+\n00000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
        let mut reader = Reader::new(seq(test));
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let r = rec.unwrap();
        assert_eq!(r.start_line_number(), 1);
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let r = rec.unwrap();
        assert_eq!(r.start_line_number(), 5);
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let r = rec.unwrap();
        assert_eq!(r.start_line_number(), 9);
    }

    #[test]
    fn test_mismatched_lengths() {
        let mut reader = Reader::new(seq(b"@test\nAGCT\n+\nIII\n@TEST\nA\n+\nI"));
        let rec = reader.next().unwrap();
        assert!(rec.is_err());
        let e = rec.unwrap_err();
        assert_eq!(e.kind, ParseErrorKind::UnequalLengths);
    }

    // https://github.com/onecodex/needletail/pull/36
    #[test]
    fn test_bad_headers() {
        let mut reader = Reader::from_path("tests/data/bad_header.fastq").unwrap();
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let rec2 = reader.next().unwrap();
        let e = rec2.unwrap_err();
        // Ideally this would not be UnexpectedEnd since we know it's an invalid record
        // but that's for another day
        assert_eq!(e.kind, ParseErrorKind::UnexpectedEnd);
    }

    // https://github.com/onecodex/needletail/pull/39
    #[test]
    fn test_fastq_with_random_tsv_inside() {
        let mut reader = Reader::from_path("tests/data/random_tsv.fq").unwrap();
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let rec2 = reader.next().unwrap();
        let e = rec2.unwrap_err();
        // It errors when it tries to validate the separator line that needs to start with `+`
        assert_eq!(e.kind, ParseErrorKind::InvalidSeparator);
    }
}
