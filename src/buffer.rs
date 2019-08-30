use std::io;
use std::marker::PhantomData;
use std::mem::transmute;

use crate::util::ParseError;
use crate::formats::fastq::{get_fastq, FASTQ};

pub struct RecBuffer {
    // rec_type: PhantomData<T>,
    file: Box<dyn io::Read>,
    pub buf: Vec<u8>,
    pub pos: usize,
    pub last: bool,
    pub count: usize,
}

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


/// A buffer that wraps an object with the `Read` trait and allows extracting
/// a set of slices to data. Acts as a lower-level primitive for our FASTX
/// readers.
impl RecBuffer {
    /// Instantiate a new buffer.
    ///
    /// # Panics
    ///
    /// Under some very rare circumstances (setting a `buf_size` larger than 2 Gb
    /// on Mac OS X) a panic can occur. Please use a smaller buffer in this case.
    pub fn new<'a>(
        mut file: &mut dyn io::Read,
        buf_size: usize,
        header: &[u8],
    ) -> Result<Self, ParseError> {
        let (buf, last) = fill_buffer(&mut file, header, buf_size)?;
        let file_ptr: *mut dyn io::Read = file;
        let file = unsafe {
            Box::from_raw(transmute(file_ptr))
        };

        Ok(RecBuffer {
            // rec_type: PhantomData,
            file,
            buf,
            pos: 0,
            last,
            count: 0,
        })
    }

    pub fn from_bytes(data: &[u8]) -> Self {
        RecBuffer {
            // rec_type: PhantomData,
            file: Box::new(io::empty()),
            buf: data.to_vec(),
            pos: 0,
            last: true,
            count: 0,
        }
    }

    /// Refill the buffer and increase its capacity if it's not big enough
    pub fn refill(&mut self) -> Result<bool, ParseError> {
        if self.pos == 0 && self.last {
            return Ok(true);
        }
        let data = &self.buf[self.pos..];
        let (buf, last) = fill_buffer(&mut self.file, &data, self.buf.capacity())?;
        self.buf = buf;
        self.last = last;
        Ok(false)
    }

    pub fn get_buffer(&self) -> Option<&[u8]> {
        if self.pos >= self.buf.len() {
            return None;
        }
        Some(&self.buf[self.pos..])
    }

    pub fn last(&self) -> bool {
        self.last
    }

    pub fn add_record(&mut self, buffer_used: usize) {
        self.pos += buffer_used;
        self.count += 1;
    }

    pub fn next<'a>(&'a mut self) -> Option<Result<FASTQ<'a>, ParseError>> {
        loop {
            if let Some(x) = get_fastq(self) {
                return Some(x);
            }
            match self.refill() {
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
