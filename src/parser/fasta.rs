//! The vast majority of the code is taken from https://github.com/markschl/seq_io/blob/master/src/fasta.rs

use crate::errors::{ErrorPosition, ParseError};
use crate::parser::record::SequenceRecord;
use crate::parser::utils::{
    fill_buf, find_line_ending, grow_to, trim_cr, FastxReader, Format, LineEnding, Position,
    BUFSIZE,
};
use memchr::{memchr2, Memchr};
use std::borrow::Cow;
use std::fs::File;
use std::{collections, io};
use std::cmp::min;
use std::io::{BufRead, Seek, SeekFrom};
use std::path::Path;

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
        println!("complete: {}", complete);

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

/// Record of a FASTA index.
#[derive(Clone, Eq, PartialEq, Debug)]
struct IndexRecord {
    name: String,
    len: u64,
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
    fetched_id: Option<IndexRecord>,
}

impl IndexedReader<File> {
    // TODO!!
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<IndexedReader<File>> {
        let reader = File::open(path).map(Reader::new).unwrap();
        let mut _inner = Vec::new();
        _inner.push(IndexRecord {
            name: String::from("xxx"),
            len: 10,
            offset: 5,
            line_bases: 5,
            line_bytes: 6,
        });
        _inner.push(IndexRecord {
            name: String::from("yyy"),
            len: 10,
            offset: 22,
            line_bases: 10,
            line_bytes: 11,
        });
        _inner.push(IndexRecord {
            name: String::from("zzz"),
            len: 2,
            offset: 38,
            line_bases: 2,
            line_bytes: 3,
        });
        let mut _name_map = collections::HashMap::new();
        _name_map.insert(String::from("xxx"), 0);
        _name_map.insert(String::from("yyy"), 1);
        _name_map.insert(String::from("zzz"), 2);
        let index = Index{ // just a dummy index
            inner: _inner,
            name_to_rid: _name_map,
        };
        Ok(IndexedReader{reader, index, start: None, end: None, fetched_id: None })
        }
}

impl<R: io::Read + Seek> IndexedReader<R> {

    pub fn fetch(&mut self, seq_name: &str, start: u64, end: u64) -> io::Result<()> {
        let rid = self.index.name_to_rid.get(seq_name).unwrap();
        let id = self.index.inner[*rid].clone();
        self.fetched_id = Some(id);
        self.start = Some(start);
        self.end = Some(end);
        Ok(())
    }

    fn seek_to(&mut self, rid: &IndexRecord, start: u64) -> io::Result<u64> {
        let crt_line_offset = start % rid.line_bases;
        println!("crt_line_offset: {}", crt_line_offset);
        let crt_line_start = start / rid.line_bases * rid.line_bytes;
        println!("crt_line_start: {}", crt_line_start);
        let crt_offset = rid.offset + crt_line_start + crt_line_offset;
        println!("crt_offset: {}", crt_offset);
        self.reader.buf_reader.seek(SeekFrom::Start(crt_offset))?;
        Ok(crt_line_offset)
    }

    fn read_into_buffer(
        &mut self,
        rid: IndexRecord,
        start: u64,
        end: u64,
        seq: &mut Vec<u8>
    ) -> io::Result<()>
    {
        let mut bases_left = end - start;
        let mut line_offset = self.seek_to(&rid, start).unwrap();
        seq.clear();
        while bases_left > 0 {
            println!("bases_left: {}", bases_left);
            let bases_read = self.read_line(&rid, &mut line_offset, bases_left, seq).unwrap();
            println!("bases_read: {}", bases_read);
            bases_left -= bases_read;
        }
        Ok(())
    }
    fn read_line(
        &mut self,
        rid: &IndexRecord,
        line_offset: &mut u64,
        bases_left: u64,
        seq: &mut Vec<u8>) -> io::Result<u64>
    {
        match fill_buf(&mut self.reader.buf_reader) {
            Ok(n) => {
                if n == 0 {
                    self.reader.finished = true;
                    // return None;
                }
            }
            Err(e) => {
                // return Some(Err(e.into()));
            }
        };
        let src = self.reader.get_buf();
        println!("src: {:?}", src);
        let rest_bases_on_crt_line = rid.line_bases - min(*line_offset, rid.line_bases);
        println!("rest_bases_on_crt_line: {}", rest_bases_on_crt_line);
        let rest_bases_in_crt_buffer = min(rest_bases_on_crt_line, src.len() as u64);
        println!("rest_bases_in_crt_buffer: {}", rest_bases_in_crt_buffer);
        let (bytes_to_read, bytes_to_keep) = if rest_bases_in_crt_buffer <= bases_left {
            let bytes_to_read = min(src.len() as u64, rid.line_bytes - *line_offset);

            (bytes_to_read, rest_bases_in_crt_buffer)
        } else {
            (bases_left, bases_left)
        };
        println!("bytes_to_read: {}, bytes_to_keep: {}", bytes_to_read, bytes_to_keep);
        seq.extend_from_slice(&src[..bytes_to_keep as usize]);

        self.reader.buf_reader.consume(bytes_to_read as usize);
        *line_offset += bytes_to_read;
        println!("line_offset: {}", line_offset);
        if *line_offset >= rid.line_bytes {
            *line_offset = 0;
        }
        Ok(bytes_to_keep)

    }

    pub fn test_fetch(&mut self, seq_name: &str, start: u64, end: u64) -> Vec<u8> {
        self.fetch(seq_name, start, end).unwrap();
        let mut seq = Vec::new();
        let rid = self.fetched_id.clone().unwrap();
        self.read_into_buffer(rid, start, end, &mut seq).unwrap();
        // println!("seq: {:?}", seq);
        seq
    }
}

/// A sub sequence of a FASTA record
#[derive(Debug)]
struct SubSequence<'a> {
    record: &'a [u8],
    name: &'a String,
    start: usize,
    end: usize,
    reverse: bool,
    complement: bool,
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
