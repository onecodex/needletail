use std::str;
//use std::io::Read;
//use std::fs::File;


//pub struct FASTXFile<'a> {
//    buffer: &'a ,
//    reader: &'a FASTXReader<'a>
//}
//
//impl<'a> FASTXFile<'a> {
//    pub fn new(file: &'a mut File) -> FASTXFile<'a> {
//        let mut buffer = Vec::new();
//        file.read_to_end(&mut buffer);
//
//        let reader = FASTXReader::new(&buffer);
//
//        FASTXFile {
//            file: file,
//            reader: &reader
//        }
//    }
//
//    pub fn reader(&'a self) -> &FASTXReader<'a> {
//        self.reader
//    }
//}

pub struct FASTXReader<'a> {
    cur_pos: usize,
    eof: bool,
    buffer: &'a [u8],
}

impl<'a> FASTXReader<'a> {
    pub fn new(slice: &'a [u8]) -> FASTXReader<'a> {
        FASTXReader {
            cur_pos: 0usize,
            eof: false,
            buffer: slice,
        }
    }

    pub fn reset_eof(&mut self) -> usize {
        //! Use this to reset the iterator (in case there was more to the buffer
        //! that we needed to read before continuing)
        // TODO: we currently return the last record in the iterator (if it had
        // anything that looked like a sequence/qual score e.g. last record)
        // so we need to backtrack and reemit that record with the rest of its
        // seq/qual before continuing (and catch the last emitted before passing
        // it along further) -> this should all be handled downstream though?
        self.eof = false;
        // also returns how much is left unread
        self.buffer.len() - self.cur_pos
    }
}

//    pub fn from_file(filename: str) {
//        let path = Path::new(filename);
//        let mut file = match File::open(&path) {
//            Err(e) => panic!("Error opening file: {}", e),
//            Ok(f) => f
//        };
//
//        FASTXReader {
//            cur_pos: Some(0usize),
//            eof: false,
//            buffer: slice
//        }
//    }
    
    // TODO: look into mmap? http://stackoverflow.com/questions/28516996/how-to-create-and-write-to-memory-mapped-files

fn find_pos(buffer: &[u8], start: usize, search: &[u8]) -> Option<usize> {
    // TODO: there has to be something like this built-in? can't find it though :/
    for (i, byte) in buffer[start..].windows(search.len()).enumerate() {
        if *byte == *search {
            return Some(i + start);
        }
    }
    None
}

impl<'a> Iterator for FASTXReader<'a> {
    type Item = (Option<&'a str>, &'a [u8], Option<&'a [u8]>);

    fn next(&mut self) -> Option<(Option<&'a str>, &'a [u8], Option<&'a [u8]>)> {
        let id_start: usize;
        match self.eof {
            true => {
                return None;
            },
            false => {
                id_start = self.cur_pos;
            }
        }
        match self.buffer[id_start] as char {
            '>' => {
                // FASTA reader
                let id_end = match find_pos(&self.buffer, id_start, &['\n' as u8]) {
                    None => {
                        self.eof = true;
                        return None
                    },
                    Some(v) => v
                };
                let seq_end = match find_pos(&self.buffer, id_end + 1, 
                                             &['\n' as u8, '>' as u8]) {
                    None => {
                        // TODO: clip off trailing newlines? theoretically anything
                        // downstream should be handling this too
                        self.eof = true;
                        self.buffer.len()
                    },
                    Some(v) => {
                        self.cur_pos = v + 1;
                        v
                    }
                };
                // clip off the '>' from the beginning of the id and make it a str
                let id = match str::from_utf8(&self.buffer[id_start + 1..id_end]) {
                    Ok(v) => Some(v),
                    Err(e) => panic!("Invalid UTF-8 in FASTA id: {}", e)
                };
                Some((id, &self.buffer[id_end+1..seq_end], None))
            },
            '@' => {
                // FASTQ reader
                let id_end = match find_pos(&self.buffer, id_start, &['\n' as u8]) {
                    None => {
                        self.eof = true;
                        return None
                    },
                    Some(v) => v
                };
                let seq_end = match find_pos(&self.buffer, id_end + 1, &['\n' as u8, '+' as u8]) {
                    None => {
                        self.eof = true;
                        return None
                    },
                    Some(v) => v
                };
                let id2_end = match find_pos(&self.buffer, seq_end + 2, &['\n' as u8]) {
                    None => {
                        self.eof = true;
                        return None
                    },
                    Some(v) => v
                };
                let qual_end = match find_pos(&self.buffer, id2_end + 1, &['\n' as u8, '@' as u8]) {
                    None => {
                        // TODO: clip off trailing newlines? theoretically anything
                        // downstream should be handling this too
                        self.eof = true;
                        self.buffer.len()
                    },
                    Some(v) => {
                        self.cur_pos = v + 1;
                        v
                    }
                };
                // clip off the '>' from the beginning of the id and make it a str
                let id = match str::from_utf8(&self.buffer[id_start + 1..id_end]) {
                    Ok(v) => Some(v),
                    Err(e) => panic!("Invalid UTF-8 in FASTA id: {}", e)
                };
                Some((id, &self.buffer[id_end + 1..seq_end], Some(&self.buffer[id2_end + 1..qual_end])))

            },
            _ => None
        }
    }
}


#[test]
fn can_parse_fastx() {
    // super minimal test examples
    let res = FASTXReader::new(b">\n").next().unwrap();
    assert!(res.0.unwrap() == "");
    assert!(res.1 == b"");

    let res = FASTXReader::new(b"@\n\n+\n").next().unwrap();
    assert!(res.0.unwrap() == "");
    assert!(res.1 == b"");
    assert!(res.2.unwrap() == b"");

    // more complex examples
    let mut fr = FASTXReader::new(b">test\nAGCT\n>test2\nGATC");
    let res = fr.next().unwrap();
    assert!(res.0.unwrap() == "test");
    assert!(res.1 == b"AGCT");
    let res = fr.next().unwrap();
    assert!(res.0.unwrap() == "test2");
    assert!(res.1 == b"GATC");

    let mut fr = FASTXReader::new(b"@test\nAGCT\n+test\nAAAA");
    let results = fr.next().unwrap();
    assert!(results.0.unwrap() == "test");
    assert!(results.1 == b"AGCT");
    assert!(results.2.unwrap() == b"AAAA");
}
