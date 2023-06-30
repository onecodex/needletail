//! The vast majority of the code is taken from https://github.com/markschl/seq_io/blob/master/src/fasta.rs

use crate::errors::{ErrorPosition, IndexError, ParseError};
use crate::parser::record::SequenceRecord;
use crate::parser::utils::{
    fill_buf, find_line_ending, grow_to, trim_cr, FastxReader, Format, LineEnding, Position,
    BUFSIZE,
};
use buffer_redux::BufReader;
use memchr::{memchr2, Memchr};
use std::borrow::Cow;
use std::cmp::min;
use std::fs::File;
use std::io::{BufRead, Seek, SeekFrom};
use std::path::Path;
use std::str;
use std::{collections, fmt, io};

#[derive(Clone, Debug)]
pub struct BufferPosition {
    /// index of '>'
    pub(crate) start: usize,
    /// Indicate line start, but actually it is one byte before (start - 1), which is usually
    /// the line terminator of the header (if there is one). The last index in the Vec is always
    /// the last byte of the last sequence line (including line terminator if present).
    /// Therefore, the length of this Vec should never be 0.
    pub(crate) seq_pos: Vec<usize>,
}

impl BufferPosition {
    #[inline]
    fn is_new(&self) -> bool {
        self.seq_pos.is_empty()
    }

    #[inline]
    fn reset(&mut self, start: usize) {
        self.seq_pos.clear();
        self.start = start;
    }

    #[inline]
    fn find_line_ending(&self, buffer: &[u8]) -> Option<LineEnding> {
        find_line_ending(self.all(buffer))
    }

    #[inline]
    pub(crate) fn all<'a>(&self, buffer: &'a [u8]) -> &'a [u8] {
        &buffer[self.start..*self.seq_pos.last().unwrap()]
    }

    #[inline]
    pub(crate) fn id<'a>(&self, buffer: &'a [u8]) -> &'a [u8] {
        trim_cr(&buffer[self.start + 1..*self.seq_pos.first().unwrap()])
    }

    #[inline]
    pub(crate) fn raw_seq<'a>(&self, buffer: &'a [u8]) -> &'a [u8] {
        if self.seq_pos.len() > 1 {
            let start = *self.seq_pos.first().unwrap() + 1;
            let end = *self.seq_pos.last().unwrap();
            trim_cr(&buffer[start..end])
        } else {
            b""
        }
    }

    #[inline]
    pub(crate) fn seq<'a>(&self, buffer: &'a [u8]) -> Cow<'a, [u8]> {
        // TODO: make that DRY
        let seq = if self.seq_pos.len() > 1 {
            let start = *self.seq_pos.first().unwrap() + 1;
            let end = *self.seq_pos.last().unwrap();
            trim_cr(&buffer[start..end])
        } else {
            b""
        };

        // first part is a fast check to see if we need to do any allocations
        let mut i;
        match memchr2(b'\r', b'\n', seq) {
            Some(break_loc) => i = break_loc,
            None => return seq.into(),
        }
        // we found a newline; create a new buffer and stripping out newlines
        // and writing into it
        let mut new_buf = Vec::with_capacity(seq.len() - 1);
        new_buf.extend_from_slice(&seq[..i]);
        while i < seq.len() {
            match memchr2(b'\r', b'\n', &seq[i..]) {
                None => {
                    new_buf.extend_from_slice(&seq[i..]);
                    break;
                }
                Some(match_pos) => {
                    new_buf.extend_from_slice(&seq[i..i + match_pos]);
                    i += match_pos + 1;
                }
            }
        }
        new_buf.into()
    }

    #[inline]
    pub(crate) fn num_bases(&self, buffer: &[u8]) -> usize {
        let seq = self.raw_seq(buffer);
        let num_lines = bytecount::count(seq, b'\n');
        let windows_num_lines = bytecount::count(seq, b'\r');
        seq.len() - num_lines - windows_num_lines
    }
}

/// Parser for FASTA files.
/// Only use this directly if you know your file is FASTA and that it is not compressed as
/// it does not handle decompression.
/// If you are unsure, it's better to use [parse_fastx_file](fn.parse_fastx_file.html).
pub struct Reader<R: io::Read> {
    buf_reader: buffer_redux::BufReader<R>,
    buf_pos: BufferPosition,
    search_pos: usize,
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
    /// use needletail::parser::{FastaReader, FastxReader};
    /// let fasta = b">id\nSEQUENCE";
    ///
    /// let mut reader = FastaReader::new(&fasta[..]);
    /// let record = reader.next().unwrap().unwrap();
    /// assert_eq!(record.id(), b"id")
    /// ```
    #[inline]
    pub fn new(reader: R) -> Reader<R> {
        Reader::with_capacity(reader, BUFSIZE)
    }

    /// Creates a new reader with a given buffer capacity. The minimum allowed
    /// capacity is 3.
    #[inline]
    pub fn with_capacity(reader: R, capacity: usize) -> Reader<R> {
        assert!(capacity >= 3);
        Reader {
            buf_reader: buffer_redux::BufReader::with_capacity(capacity, reader),
            buf_pos: BufferPosition {
                start: 0,
                seq_pos: Vec::with_capacity(1),
            },
            position: Position::new(0, 0),
            search_pos: 0,
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
    /// use needletail::parser::{FastaReader, FastxReader};
    ///
    /// let mut reader = FastaReader::from_path("seqs.fasta").unwrap();
    ///
    /// // (... do something with the reader)
    /// ```
    #[inline]
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Reader<File>> {
        File::open(path).map(Reader::new)
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

    #[inline]
    fn next_pos(&mut self) {
        self.position.line += self.buf_pos.seq_pos.len() as u64;
        self.position.byte += (self.search_pos - self.buf_pos.start) as u64;
        self.buf_pos.reset(self.search_pos);
    }

    /// Finds the position of the next record
    /// and returns true if found; false if end of buffer reached.
    #[inline]
    fn find(&mut self) -> bool {
        if self._find() {
            return true;
        }

        // nothing found
        if self.get_buf().len() < self.buf_reader.capacity() {
            // EOF reached, there will be no next record
            self.finished = true;
            if !self.buf_pos.seq_pos.is_empty() {
                self.buf_pos.seq_pos.push(self.search_pos);
            }
            return true;
        }

        false
    }

    /// Returns true if complete position found, false if end of buffer reached.
    #[inline]
    fn _find(&mut self) -> bool {
        let bufsize = self.get_buf().len();
        for pos in Memchr::new(b'\n', &self.buf_reader.buffer()[self.search_pos..]) {
            let pos = self.search_pos + pos;
            let next_line_start = pos + 1;

            if next_line_start == bufsize {
                // cannot check next byte -> treat as incomplete
                self.search_pos = pos; // make sure last byte is re-searched next time
                return false;
            }

            self.buf_pos.seq_pos.push(pos);
            if self.get_buf()[next_line_start] == b'>' {
                // complete record was found
                self.search_pos = next_line_start;
                return true;
            }
        }

        // record end not found
        self.search_pos = bufsize;
        false
    }

    /// To be called when the end of the buffer is reached and `next_pos` does not find
    /// the next record. Incomplete bytes will be moved to the start of the buffer.
    /// If the record still doesn't fit in, the buffer will be enlarged.
    /// After calling this function, the position will therefore always be 'complete'.
    /// this function assumes that the buffer was fully searched
    fn next_complete(&mut self) -> Result<bool, ParseError> {
        loop {
            if self.buf_pos.start == 0 {
                // first record -> buffer too small
                self.grow();
            } else {
                // not the first record -> buffer may be big enough
                self.make_room();
            }

            // fill up remaining buffer
            fill_buf(&mut self.buf_reader)?;

            if self.find() {
                return Ok(true);
            }
        }
    }

    /// Grow internal buffer as needed
    fn grow(&mut self) {
        let cap = self.buf_reader.capacity();
        let new_size = grow_to(cap);
        let additional = new_size - cap;
        self.buf_reader.reserve(additional);
    }

    /// Move incomplete bytes to start of buffer
    fn make_room(&mut self) {
        let consumed = self.buf_pos.start;
        self.buf_reader.consume(consumed);
        self.buf_reader.make_room();
        self.buf_pos.start = 0;
        self.search_pos -= consumed;
        for s in &mut self.buf_pos.seq_pos {
            *s -= consumed;
        }
    }
}

impl<R: io::Read + Send> FastxReader for Reader<R> {
    fn next(&mut self) -> Option<Result<SequenceRecord, ParseError>> {
        if self.finished {
            return None;
        }

        // Load some data in the buffer to start
        if self.position.line == 0 {
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

            if self.get_buf()[0] == b'>' {
                self.position.line = 1;
                self.position.byte = 0;
                self.buf_pos.start = 0;
                self.search_pos = 1;
            } else {
                return Some(Err(ParseError::new_invalid_start(
                    self.get_buf()[0],
                    ErrorPosition {
                        line: self.position.line,
                        id: None,
                    },
                    Format::Fasta,
                )));
            }
        }

        if !self.buf_pos.is_new() {
            self.next_pos();
        }

        // Can we identify the start of the next record ?
        let complete = self.find();

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

        if self.buf_pos.seq_pos.is_empty() {
            return Some(Err(ParseError::new_unexpected_end(
                ErrorPosition {
                    line: self.position.line,
                    id: None,
                },
                Format::Fasta,
            )));
        }

        if self.line_ending.is_none() {
            self.line_ending = self.buf_pos.find_line_ending(self.get_buf());
        }
        Some(Ok(SequenceRecord::new_fasta(
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

/// A FASTA index including info at..
#[derive(Debug)]
struct Index {
    inner: Vec<IndexRecord>,
    name_to_rid: collections::HashMap<String, usize>,
}

impl Index {
    /// Create a new index from a faidx file.
    fn new<R: io::Read>(fai: R) -> Result<Index, IndexError> {
        let mut inner = Vec::new();
        let mut name_to_rid = collections::HashMap::new();
        for (rid, line) in BufReader::new(fai).lines().flatten().enumerate() {
            let values: Vec<&str> = line.split(char::is_whitespace).collect();
            if values.len() != 5 {
                return Err(IndexError::new_fai_format_err());
            }
            let name = values[0].to_owned();
            let indexrecord = IndexRecord {
                name,
                len: values[1].parse::<u64>()?, // convert intErr to ParseError
                offset: values[2].parse::<u64>()?,
                line_bases: values[3].parse::<u64>()?,
                line_bytes: values[4].parse::<u64>()?,
            };
            name_to_rid.insert(indexrecord.name.clone(), rid);
            inner.push(indexrecord);
        }
        Ok(Index { inner, name_to_rid })
    }

    /// Open a FASTA index file from a given path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Result<Index, IndexError>> {
        File::open(path).map(Index::new)
    }

    /// Infer the path to the index file from the path to the given FASTA file.
    pub fn from_fasta_path<P: AsRef<Path>>(path: P) -> io::Result<Result<Index, IndexError>> {
        let mut path = path.as_ref().as_os_str().to_owned();
        path.push(".fai");
        Index::from_path(path)
    }
}

/// Record of a FASTA index.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct IndexRecord {
    pub name: String,
    pub len: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

/// A FASTER reader with a index
pub struct IndexedReader<R: io::Read> {
    reader: Reader<R>,
    index: Index,
    start: Option<u64>,
    end: Option<u64>,
    pub fetched_id: Option<IndexRecord>,
}

impl IndexedReader<File> {
    /// Read a FASTA file and its index from a given path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<IndexedReader<File>, IndexError> {
        // Judge fai file exist firstly
        let index = match Index::from_fasta_path(path.as_ref()) {
            Ok(Ok(index)) => index,
            Ok(Err(e)) => return Err(e),
            _ => return Err(IndexError::new_fai_io_err(path.as_ref().display())),
        };
        let reader = Reader::from_path(path)?;
        Ok(IndexedReader {
            reader,
            index,
            start: None,
            end: None,
            fetched_id: None,
        })
    }
}

impl<R: io::Read + Seek> IndexedReader<R> {
    /// fetch to a given region
    pub fn fetch(
        &mut self,
        seq_name: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> Result<(), IndexError> {
        let rid = match self.index.name_to_rid.get(seq_name) {
            Some(rid) => rid,
            None => return Err(IndexError::new_seq_name_err(seq_name)),
        };
        let id = self.index.inner[*rid].clone(); // rid should be in inner
        self.fetched_id = Some(id);
        self.start = start;
        self.end = end;
        Ok(())
    }

    /// seek to the position of the given record and start
    fn seek_to(&mut self, rid: &IndexRecord, start: u64) -> io::Result<u64> {
        let crt_line_offset = start % rid.line_bases; // current line offset
        let crt_line_start = start / rid.line_bases * rid.line_bytes; // current line start
        let crt_offset = rid.offset + crt_line_start + crt_line_offset; // current offset
        self.reader.buf_reader.seek(SeekFrom::Start(crt_offset))?;
        Ok(crt_line_offset)
    }

    /// read sub sequence into a buffer
    fn read_into_buffer(
        &mut self,
        rid: &IndexRecord,
        start: u64,
        end: u64,
        seq: &mut Vec<u8>,
    ) -> Result<(), IndexError> {
        // initialize
        let mut bases_rest = end - start;
        let mut line_offset = self.seek_to(rid, start)?;
        seq.clear();
        // start reading
        while bases_rest > 0 {
            let bases_read = self.read_line(rid, &mut line_offset, bases_rest, seq)?;
            bases_rest -= bases_read;
        }
        Ok(())
    }

    /// read a line into a buffer
    fn read_line(
        &mut self,
        rid: &IndexRecord,
        line_offset: &mut u64,
        bases_rest: u64,
        seq: &mut Vec<u8>,
    ) -> Result<u64, IndexError> {
        // fill buffer and get buffer
        let src = self.reader.buf_reader.fill_buf()?;

        // get rest bases on current line
        let rest_bases_on_crt_line = rid.line_bases - min(*line_offset, rid.line_bases);
        // get rest bases on current buffer
        let rest_bases_on_crt_buffer = min(rest_bases_on_crt_line, src.len() as u64);

        // compute bytes to read and bytes to keep
        let (bytes_to_read, bytes_to_keep) = if rest_bases_on_crt_buffer <= bases_rest {
            let bytes_to_read = min(src.len() as u64, rid.line_bytes - *line_offset);

            (bytes_to_read, rest_bases_on_crt_buffer)
        } else {
            (bases_rest, bases_rest)
        };

        // extend seq from src slice
        seq.extend_from_slice(&src[..bytes_to_keep as usize]);

        // consume bytes_to_read in buffer
        self.reader.buf_reader.consume(bytes_to_read as usize);

        // move line_offset
        *line_offset += bytes_to_read;
        // stop at the end of line
        if *line_offset >= rid.line_bytes {
            *line_offset = 0;
        }

        Ok(bytes_to_keep)
    }

    pub fn subseq(
        &mut self,
        seq_name: &str,
        start: Option<u64>,
        end: Option<u64>,
    ) -> Result<SubSequence, IndexError> {
        // fetch to the given region
        self.fetch(seq_name, start, end)?;

        // get the rid and the true start and end
        let rid = self.fetched_id.clone().unwrap(); // a safe unwrap
        let start = self.start.unwrap_or(0); // None -> 0
        let end = self.end.unwrap_or(rid.len); // None -> length

        // check if the region is valid
        if start > end || end > rid.len {
            // don't judge if start < 0 due to u64
            return Err(IndexError::new_invalid_region_err());
        }

        // prepare the buffer
        let mut seq = Vec::new();

        // read into the buffer
        self.read_into_buffer(&rid, start, end, &mut seq)?;

        // return the subsequence
        Ok(SubSequence {
            seq,
            name: rid.name,
            start,
            end,
        })
    }
}

/// A sub sequence of a FASTA record
pub struct SubSequence {
    pub seq: Vec<u8>,
    pub name: String,
    pub start: u64,
    pub end: u64,
}

impl fmt::Display for SubSequence {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let seq_str = str::from_utf8(&self.seq).unwrap();
        write!(f, ">{}:{}-{}\n{}", self.name, self.start, self.end, seq_str)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use super::*;
    use crate::errors::ParseErrorKind;

    fn seq(s: &[u8]) -> Cursor<&[u8]> {
        Cursor::new(&s[..])
    }

    #[test]
    fn test_basic() {
        let mut reader = Reader::new(seq(b">test\nACGT\n>test2\nTGCA\n"));
        assert!(reader.line_ending().is_none());
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let r = rec.unwrap();
        assert_eq!(r.id(), b"test");
        assert_eq!(r.raw_seq(), b"ACGT");
        assert_eq!(r.all(), b">test\nACGT");
        assert_eq!(reader.line_ending().unwrap(), LineEnding::Unix);
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let r = rec.unwrap();
        assert_eq!(r.id(), b"test2");
        assert_eq!(r.raw_seq(), b"TGCA");
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_wrapped_fasta() {
        let mut reader = Reader::new(seq(b">test\nACGT\nACGT\n>test2\nTGCA\nTG"));
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let r = rec.unwrap();
        assert_eq!(r.id(), b"test");
        assert_eq!(r.raw_seq(), b"ACGT\nACGT");
        assert_eq!(r.num_bases(), 8);
        assert_eq!(reader.line_ending().unwrap(), LineEnding::Unix);
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let r = rec.unwrap();
        assert_eq!(r.id(), b"test2");
        assert_eq!(r.raw_seq(), b"TGCA\nTG");
        assert_eq!(r.num_bases(), 6);
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_wrapped_fasta_windows_newlines() {
        let mut reader = Reader::new(seq(b">test\r\nACGT\r\nACGT\r\n>test2\r\nTGCA\r\nTG"));
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let r = rec.unwrap();
        assert_eq!(r.id(), b"test");
        assert_eq!(r.raw_seq(), b"ACGT\r\nACGT");
        assert_eq!(r.num_bases(), 8);
        assert_eq!(r.start_line_number(), 1);
        assert_eq!(reader.line_ending().unwrap(), LineEnding::Windows);
        let rec = reader.next().unwrap();
        assert!(rec.is_ok());
        let r = rec.unwrap();
        assert_eq!(r.id(), b"test2");
        assert_eq!(r.raw_seq(), b"TGCA\r\nTG");
        assert_eq!(r.num_bases(), 6);
        assert_eq!(r.start_line_number(), 4);
        assert!(reader.next().is_none());
    }

    #[test]
    fn test_premature_ending() {
        let mut reader = Reader::new(seq(b">test\nAGCT\n>test2"));
        reader.next().unwrap().unwrap();
        let rec = reader.next().unwrap();
        assert!(rec.is_err());
        let r = rec.unwrap_err();
        assert_eq!(r.kind, ParseErrorKind::UnexpectedEnd);

        let mut reader = Reader::new(seq(b">test\r\nAGCT\r\n>test2\r\n"));
        reader.next().unwrap().unwrap();
        let rec = reader.next().unwrap();
        assert!(rec.is_err());
        let r = rec.unwrap_err();
        assert_eq!(r.kind, ParseErrorKind::UnexpectedEnd);
    }

    #[test]
    fn test_empty_records() {
        let mut reader = Reader::new(seq(b">\n\n>shine\nAGGAGGU"));
        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.id(), b"");
        assert_eq!(rec.raw_seq(), b"");
        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.id(), b"shine");
        assert_eq!(rec.raw_seq(), b"AGGAGGU");

        let mut reader = Reader::new(seq(b">\r\n\r\n>shine\r\nAGGAGGU"));
        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.id(), b"");
        assert_eq!(rec.raw_seq(), b"");
        let rec = reader.next().unwrap().unwrap();
        assert_eq!(rec.id(), b"shine");
        assert_eq!(rec.raw_seq(), b"AGGAGGU");
    }
}
