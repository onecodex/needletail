use std::io;
use std::marker::PhantomData;

use crate::util::ParseError;

#[inline]
fn fill_buffer(
    file: &mut dyn io::Read,
    data: &[u8],
    buf_size: usize,
) -> Result<(Vec<u8>, bool), ParseError> {
    let mut buf = Vec::with_capacity(buf_size + data.len());
    unsafe {
        buf.set_len(buf_size + data.len());
    }
    buf[..data.len()].copy_from_slice(data);
    let amt_read = file.read(&mut buf[data.len()..])?;
    unsafe {
        buf.set_len(amt_read + data.len());
    }
    Ok((buf, amt_read == 0))
}

pub struct RecBuffer<'a, T> {
    rec_type: PhantomData<T>,
    file: Option<&'a mut dyn io::Read>,
    pub buf: Vec<u8>,
    pub last: bool,
}

/// A buffer that wraps an object with the `Read` trait and allows extracting
/// a set of slices to data. Acts as a lower-level primitive for our FASTX
/// readers.
impl<'a, 'b, T> RecBuffer<'a, T>
where
    T: RecReader<'b>,
{
    /// Instantiate a new buffer.
    ///
    /// # Panics
    ///
    /// Under some very rare circumstances (setting a `buf_size` larger than 2 Gb
    /// on Mac OS X) a panic can occur. Please use a smaller buffer in this case.
    pub fn new(
        file: &'a mut dyn io::Read,
        buf_size: usize,
        header: &[u8],
    ) -> Result<RecBuffer<'a, T>, ParseError> {
        let (buf, last) = fill_buffer(file, header, buf_size)?;
        Ok(RecBuffer {
            rec_type: PhantomData,
            file: Some(file),
            last,
            buf,
        })
    }

    pub fn new_chunked() -> Result<RecBuffer<'a, T>, ParseError> {
        Ok(RecBuffer {
            rec_type: PhantomData,
            file: None,
            last: false,
            buf: Vec::new(),
        })
    }

    pub fn fill(&mut self, data: &[u8], last: bool) -> Result<(), ParseError> {
        let mut data = io::Cursor::new(data);
        let (buf, _) = fill_buffer(&mut data, &self.buf, self.buf.capacity())?;
        self.buf = buf;
        self.last = last;
        Ok(())
    }

    /// Refill the buffer and increase its capacity if it's not big enough.
    /// Takes a tuple of the bytes used and how many records returned so far.
    #[inline]
    pub fn refill(&mut self, used: usize) -> Result<bool, ParseError> {
        if used == 0 && self.last {
            return Ok(true);
        }
        let data = &self.buf[used..];
        let (buf, last) = if let Some(f) = &mut self.file {
            fill_buffer(f, &data, self.buf.capacity())?
        } else {
            (data.to_vec(), self.last)
        };
        self.buf = buf;
        self.last = last;
        Ok(false)
    }

    pub fn get_reader(&'b self) -> T {
        T::from_buffer(&self.buf, self.last)
    }
}

/// RecReader is an adaptor trait that allows new file format parsers to be
/// defined. It takes a chunk from a RecBuffer (`from_reader`), optionally
/// parses an initial header out (`header`) and then provides an iterator
/// interface to parse a record stream. When finished, it provides a `eof`
/// function to determine if the stream is completely exhausted.
///
pub trait RecReader<'s>: Sized + Iterator {
    type Header;

    fn from_buffer(buf: &'s [u8], last: bool) -> Self;
    fn header(&mut self) -> Result<Self::Header, ParseError>;
    fn eof(&self) -> Result<(), ParseError>;
    fn used(&self) -> usize;
}
