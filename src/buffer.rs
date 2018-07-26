use std::io;
use std::marker::PhantomData;

use util::ParseError;


pub struct RecReader<'a> {
    file: &'a mut io::Read,
    last: bool,
    buf: Vec<u8>,
}

/// A buffer that wraps an object with the `Read` trait and allows extracting
/// a set of slices to data. Acts as a lower-level primitive for our FASTX
/// readers.
impl<'a> RecReader<'a> {
    /// Instantiate a new buffer.
    ///
    /// # Panics
    ///
    /// Under some very rare circumstances (setting a `buf_size` larger than 2 Gb
    /// on Mac OS X) a panic can occur. Please use a smaller buffer in this case.
    pub fn new(file: &'a mut io::Read, buf_size: usize, header: &[u8]) -> Result<RecReader<'a>, ParseError> {
        let mut buf = Vec::with_capacity(buf_size);
        unsafe {
            buf.set_len(buf_size + header.len());
        }
        buf[..header.len()].copy_from_slice(header);
        let amt_read = file.read(&mut buf[header.len()..])?;
        unsafe {
            buf.set_len(amt_read + header.len());
        }
        
        Ok(RecReader {
            file: file,
            last: false,
            buf: buf,
        })
    }

    /// Refill the buffer and increase its capacity if it's not big enough
    pub fn refill(&mut self, used: usize) -> Result<bool, ParseError> {
        if used == 0 && self.last {
            return Ok(true);
        }
        // if used >= self.buf.len() {
        //     return Ok(true);
        // }
        let cur_length = self.buf.len() - used;
        let new_length = cur_length + self.buf.capacity();

        let mut new_buf = Vec::with_capacity(new_length);
        unsafe {
            new_buf.set_len(new_length);
        }
        new_buf[..cur_length].copy_from_slice(&self.buf[used..]);
        let amt_read = self.file.read(&mut new_buf[cur_length..])?;
        unsafe {
            new_buf.set_len(cur_length + amt_read);
        }
        self.buf = new_buf;
        self.last = amt_read == 0;
        Ok(false)
    }

    pub fn get_buffer<'b, T>(&'b self) -> RecBuffer<'b, T> {
        RecBuffer {
            buf: &self.buf,
            pos: 0,
            last: self.last,
            record_type: PhantomData,
        }
    }
}

pub struct RecBuffer<'a, T> {
    pub buf: &'a [u8],
    pub pos: usize,
    pub last: bool,
    record_type: PhantomData<T>,
}

impl<'a, T> RecBuffer<'a, T> {
    pub fn from_bytes(data: &'a [u8]) -> Self {
        RecBuffer {
            buf: data,
            pos: 0,
            last: false,
            record_type: PhantomData,
        }
    }
}

pub trait FindRecord {
    fn move_to_next(&mut self);
    fn is_finished(&self) -> bool;
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
