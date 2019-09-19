use std::io;

use safemem::copy_over;

use crate::util::ParseError;

/// A buffer that wraps an object with the `Read` trait and allows extracting
/// a set of slices to data. Acts as a lower-level primitive for our FASTX
/// readers.
pub struct RecBuffer<'a> {
    file: &'a mut dyn io::Read,
    pub buf: Vec<u8>,
    pub last: bool,
}

impl<'a> RecBuffer<'a> {
    /// Instantiate a new buffer.
    ///
    /// # Panics
    ///
    /// Under some very rare circumstances (setting a `buf_size` larger than 2 Gb
    /// on Mac OS X) a panic can occur. Please use a smaller buffer in this case.
    pub fn new(file: &'a mut dyn io::Read, buf: Vec<u8>) -> Result<RecBuffer<'a>, ParseError> {
        Ok(RecBuffer {
            file,
            last: false,
            buf,
        })
    }

    /// Refill the buffer and increase its capacity if it's not big enough.
    /// Takes a tuple of the bytes used and how many records returned so far.
    #[inline]
    pub fn refill(&mut self, used: usize) -> Result<bool, ParseError> {
        if used == 0 && self.last {
            return Ok(true);
        }
        let remaining = self.buf.len() - used;
        if used == 0 {
            let mut new_buf = Vec::with_capacity(2 * self.buf.len());
            unsafe {
                new_buf.set_len(new_buf.capacity());
            }
            new_buf[..self.buf.len()].copy_from_slice(&self.buf);
            self.buf = new_buf;
        } else if remaining != 0 {
            copy_over(&mut self.buf, used, 0, remaining);
        }
        let amt_read = self.file.read(&mut self.buf[remaining..])?;
        unsafe {
            self.buf.set_len(remaining + amt_read);
        }
        self.last = amt_read == 0;
        Ok(false)
    }
}

/// [⚠️Unstable] RecParser is an adaptor trait that allows new file format
/// parsers to be defined. It takes a chunk from a RecBuffer (`from_reader`),
/// optionally parses an initial header out (`header`) and then provides an
/// iterator interface to parse a record stream. When finished, it provides an
/// `eof` function to determine if the stream is completely exhausted.
pub trait RecParser<'s>: Sized + Iterator {
    type Header;

    fn from_buffer(buf: &'s [u8], last: bool) -> Self;
    fn header(&mut self) -> Result<Self::Header, ParseError>;
    fn eof(&self) -> Result<(), ParseError>;
    fn used(&self) -> usize;
}
