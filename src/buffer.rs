use std::io;

use crate::util::ParseError;

pub struct RecBuffer<'a> {
    file: Option<&'a mut dyn io::Read>,
    pub buf: Vec<u8>,
    pub pos: usize,
    pub last: bool,
    pub count: usize,
}

/// A buffer that wraps an object with the `Read` trait and allows extracting
/// a set of slices to data. Acts as a lower-level primitive for our FASTX
/// readers.
impl<'a> RecBuffer<'a> {
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
    ) -> Result<Self, ParseError> {
        let mut buf = Vec::with_capacity(buf_size + header.len());
        unsafe {
            buf.set_len(buf_size + header.len());
        }
        buf[..header.len()].copy_from_slice(header);
        let amt_read = file.read(&mut buf[header.len()..])?;
        unsafe {
            buf.set_len(amt_read + header.len());
        }

        Ok(RecBuffer {
            file: Some(file),
            buf,
            pos: 0,
            last: false,
            count: 0,
        })
    }

    pub fn from_bytes(data: &'a [u8]) -> Self {
        RecBuffer {
            file: None,
            buf: data.to_vec(),
            pos: 0,
            last: true,
            count: 0,
        }
    }

    /// Refill the buffer and increase its capacity if it's not big enough
    pub fn refill(&mut self, data: &[u8]) -> Result<bool, ParseError> {
        if self.pos == 0 && self.last {
            return Ok(true);
        }
        let cur_length = self.buf.len() - self.pos;
        let new_length = cur_length + self.buf.capacity();

        let mut new_buf = Vec::with_capacity(new_length);
        unsafe {
            new_buf.set_len(new_length);
        }
        new_buf[..cur_length].copy_from_slice(&self.buf[self.pos..]);
        new_buf[cur_length..cur_length + data.len()].copy_from_slice(data);
        let amt_read = if let Some(file) = &mut self.file {
            file.read(&mut new_buf[cur_length + data.len()..])?
        } else {
            0
        };
        unsafe {
            new_buf.set_len(cur_length + data.len() + amt_read);
        }
        self.buf = new_buf;
        self.last = amt_read == 0;
        Ok(false)
    }

    pub fn next<'s: 'b, 'b, T>(&'s mut self) -> Option<Result<T, ParseError>>
    where
        T: RecordFormat<'b>,
    {
        loop {
            if let Some(x) = T::parse(self) {
                return Some(x);
            }
            match self.refill(b"") {
                Err(e) => return Some(Err(e)),
                Ok(true) => return None,
                Ok(false) => {},
            };
        }
    }
}

#[test]
fn test_from_bytes() {
    // this is not a great test
    let rb: RecBuffer = RecBuffer::from_bytes(b"test");
    assert_eq!(rb.pos, 0);
    assert_eq!(rb.buf, b"test");
}

pub trait RecordFormat<'b> {
    fn parse(rbuf: &'b mut RecBuffer) -> Option<Result<Self, ParseError>>
    where
        Self: Sized;
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
