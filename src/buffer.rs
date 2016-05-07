use std::error;
use std::fmt;
use std::io;


/// Errors returned during parsing/reading buffers
#[derive(Clone, Debug, PartialEq)]
pub enum ParseError {
    NeedMore,
    EOF,
    Invalid(String),
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Parse Error: {}", *self)
    }
}

impl error::Error for ParseError {
    fn description(&self) -> &str {
        "ParseError"
    }

    fn cause(&self) -> Option<&error::Error> {
        None
    }
}


impl From<io::Error> for ParseError {
    fn from(err: io::Error) -> ParseError {
        ParseError::Invalid(err.to_string())
    }
}


pub struct RecBuffer<'a> {
    file: &'a mut io::Read,
    buf: Vec<u8>,
    offset: usize,
    record_start: usize,
    marks: Vec<usize>,
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
    pub fn new(file: &'a mut io::Read, buf_size: usize) -> RecBuffer {
        let mut buf = Vec::with_capacity(buf_size);
        unsafe {
            buf.set_len(buf_size);
        }
        let amt_read = file.read(&mut buf).unwrap();
        unsafe {
            buf.set_len(amt_read);
        }
        
        RecBuffer {
            file: file,
            buf: buf,
            offset: 0,
            record_start: 0,
            marks: Vec::new(),
        }
    }

    /// Internal method to refill the internal buffer (and optionally increase
    /// it's capacity if it's not big enough)
    fn refresh(&mut self) -> Result<bool, ParseError> {
        let cur_length = self.buf.len() - self.record_start;
        let new_length = cur_length + self.buf.capacity();

        let mut new_buf = Vec::with_capacity(new_length);
        unsafe {
            new_buf.set_len(new_length);
        }
        new_buf[..cur_length].copy_from_slice(&self.buf[self.record_start..]);
        let amt_read = self.file.read(&mut new_buf[cur_length..])?;
        unsafe {
            new_buf.set_len(cur_length + amt_read);
        }
        self.buf = new_buf;
        self.offset -= self.record_start;
        self.record_start = 0;
        
        // return eof if we didn't get anything
        Ok(amt_read == 0)
    }

    pub fn mark_field<F>(&mut self, rec_fn: F) -> Result<usize, ParseError>
        where F: Fn(&[u8], bool) -> Result<usize, ParseError>
    {
        let mut eof = false;
        loop {
            if self.offset == self.buf.len() {
                eof = self.refresh()?;
                if eof {
                    return Err(ParseError::EOF);
                }
            }
            let result = rec_fn(&self.buf[self.offset..], eof);
            match (result, eof) {
                (Ok(v), _) => {
                    self.marks.push(v);
                    self.offset += v;
                    return Ok(v);
                },
                (Err(ParseError::NeedMore), false) => {
                    eof = self.refresh()?;
                },
                (Err(e), false) => return Err(e),
                _ => return Err(ParseError::EOF),
            }
        }
    }

    pub fn fields(&mut self) -> Vec<&[u8]> {
        let mut slices = Vec::with_capacity(self.marks.len());
        let mut cum_pos = self.record_start;
        for mark in &self.marks {
            slices.push(&self.buf[cum_pos..cum_pos + *mark]);
            cum_pos += *mark;
        }

        self.marks.clear();
        self.record_start = self.offset;
        slices
    }
}



#[test]
fn test_buffer() {
    let test: [u8; 4] = [1, 2, 3, 4];
    let mut slice = &test[..];
    let mut b = RecBuffer::new(&mut slice, 2usize);

    let field_len = b.mark_field(|buf, eof| {
        assert_eq!(buf[0], 1);
        assert_eq!(buf.len(), 2);
        assert!(!eof);
        Ok(1)
    });
    assert_eq!(field_len.unwrap(), 1);

    let field_len = b.mark_field(|buf, eof| {
        if buf.len() == 1 {
            assert_eq!(buf[0], 2);
            assert!(!eof);
            return Err(ParseError::NeedMore);
        }
        if !eof {
            assert_eq!(buf[0], 2);
            assert_eq!(buf.len(), 3);
            return Err(ParseError::NeedMore);
        }

        assert_eq!(buf[0], 2);
        assert!(eof);
        assert_eq!(buf.len(), 3);
        Ok(3)
    });
    assert_eq!(field_len.unwrap(), 3);

    let records = b.fields();
    assert_eq!(records[0], &[1]);
    assert_eq!(records[1], &[2, 3, 4]);
}
