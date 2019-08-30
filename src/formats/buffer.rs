use std::io;
use std::marker::PhantomData;

use crate::util::ParseError;

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

pub struct RecReader<'a, T> {
    rec_type: PhantomData<T>,
    file: Option<&'a mut dyn io::Read>,
    pub buf: Vec<u8>,
    last: bool,
    count: usize,
}

/// A buffer that wraps an object with the `Read` trait and allows extracting
/// a set of slices to data. Acts as a lower-level primitive for our FASTX
/// readers.
impl<'a, T> RecReader<'a, T> {
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
    ) -> Result<RecReader<'a, T>, ParseError> {
        let (buf, last) = fill_buffer(file, header, buf_size)?;
        Ok(RecReader {
            rec_type: PhantomData,
            file: Some(file),
            last,
            buf,
            count: 0,
        })
    }

    pub fn new_chunked() -> Result<RecReader<'a, T>, ParseError> {
        Ok(RecReader {
            rec_type: PhantomData,
            file: None,
            last: false,
            buf: Vec::new(),
            count: 0,
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
    pub fn refill(&mut self, used: (usize, usize)) -> Result<bool, ParseError> {
        if used.0 == 0 && self.last {
            return Ok(true);
        }
        let data = &self.buf[used.0..];
        let (buf, last) = if let Some(f) = &mut self.file {
            fill_buffer(f, &data, self.buf.capacity())?
        } else {
            (data.to_vec(), self.last)
        };
        self.buf = buf;
        self.last = last;
        self.count += used.1;
        Ok(false)
    }

    pub fn get_buffer<'b>(&'b self) -> RecBuffer<'b, T> {
        RecBuffer {
            buf: &self.buf,
            pos: 0,
            last: self.last,
            record_type: PhantomData,
            count: self.count,
        }
    }
}

#[derive(Debug)]
pub struct RecBuffer<'a, T> {
    record_type: PhantomData<T>,
    pub buf: &'a [u8],
    pub pos: usize,
    pub last: bool,
    pub count: usize,
}

impl<'a, T> RecBuffer<'a, T> {
    pub fn from_bytes(data: &'a [u8]) -> Self {
        RecBuffer {
            buf: data,
            pos: 0,
            last: true,
            record_type: PhantomData,
            count: 0,
        }
    }

    pub fn used(&self) -> (usize, usize) {
        (self.pos, self.count)
    }
}

#[test]
fn test_from_bytes() {
    // this is not a useful test, but it does get the compiler to shut up
    // about `from_bytes` not being used
    let rb: RecBuffer<String> = RecBuffer::from_bytes(b"test");
    assert_eq!(rb.pos, 0);
    assert_eq!(rb.buf, b"test");
}

// pub fn parse<T, E, F>(reader: &'s mut io::Read, header: &[u8], ref mut callback: F) -> Result<(), E> where
//     E: From<ParseError>,
//     F: FnMut(T) -> Result<(), E>,
//     for<'s> RecBuffer<'s, T>: Iterator<Item=Result<T, ParseError>>,
// {
//     let mut rec_reader = RecReader::new(reader, 10_000_000, header)?;
//     loop {
//         let used = {
//             let mut rec_buffer = rec_reader.get_buffer();
//             for s in rec_buffer.by_ref() {
//                 callback(s?)?;
//             }
//             rec_buffer.pos
//         };
//         if rec_reader.refill(used)? {
//             break;
//         }
//     }
//     if rec_reader.get_buffer::<T>().is_finished(true) {
//         Ok(())
//     } else {
//         Err(ParseError::PrematureEOF.into())
//     }
// }
