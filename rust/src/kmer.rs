use std::borrow::Cow;
// TODO: should have IKmers

pub fn complement(n: &u8) -> u8 {
    match *n as char {
        'a' => 't' as u8,
        't' => 'a' as u8,
        'g' => 'c' as u8,
        'c' => 'g' as u8,
        'A' => 'T' as u8,
        'T' => 'A' as u8,
        'G' => 'C' as u8,
        'C' => 'G' as u8,
        // TODO: actually implement logic for IUPAC bases here?
        x => x as u8,
    }
}

pub fn canonical<'a>(seq: Cow<'a, [u8]>) -> Cow<'a, [u8]> {
    let mut buf: Vec<u8> = Vec::with_capacity(seq.len());
    // enough just keeps our comparisons from happening after they need to
    let mut enough = false;
    let mut original_was_canonical = false;

    // loop through the kmer and its reverse complement simultaneously
    for (rn, n) in seq.iter().rev().map(|n| complement(n)).zip(seq.iter()) {
        buf.push(rn);
        if !enough && n < &rn {
            original_was_canonical = true;
            break;
        } else if !enough && &rn < n {
            enough = true;
        }
        // unstated if branch: if rn == n, keep comparing
    }
    match (original_was_canonical, enough) {
        (true, true) => panic!("Bug: should never set original_was_canonical if enough == true"),
        (true, false) => seq,
        (false, true) => Cow::Owned(buf),
        (false, false) => seq  // the sequences were completely equal, return the ref
    }
}

#[test]
fn can_canonicalize() {
    // TODO: figure out a way to compare a Cow::Owned?
    //assert!(canonical(Cow::Borrowed(b"A")) == Cow::Owned::<[u8]>(b"T".to_vec()));
    assert!(canonical(Cow::Borrowed(b"AAGT")) == Cow::Borrowed(b"AAGT"));
    assert!(canonical(Cow::Borrowed(b"GC")) == Cow::Borrowed(b"GC"));
}

// TODO
//pub fn skip_n<'a, T>(iter: T) -> T where T: Iterator<Item=&'a [u8]> {
//    iter.filter(|kmer| kmer.contains(&('N' as u8)) || kmer.contains(&('n' as u8)))
//}

fn has_no_n(seq: Cow<[u8]>) -> bool {
    seq.iter().all(|n| match *n as char {
        'a' | 'c' | 'g' | 't' | 'A' | 'C' | 'G' | 'T' => true,
        _ => false
    })
}

pub struct Kmer<'a> {
    k: usize,
    multiline: bool,
    cur_pos: Option<usize>,
    buffer: &'a [u8],
}

impl<'a> Kmer<'a> {
    pub fn new(slice: &'a [u8], k: usize, multiline: bool) -> Kmer<'a> {
        Kmer {
            k: k,
            multiline: multiline,
            cur_pos: Some(0),
            buffer: slice,
        }
    }
}

impl<'a> Iterator for Kmer<'a> {
    type Item = Cow<'a, [u8]>;

    fn next(&mut self) -> Option<Cow<'a, [u8]>> {
        match self.cur_pos {
                    
            // the sequence is exhausted, return None as the Iterator sentinel value
            None => None,
            Some(mut pos) => {
                if pos > self.buffer.len() - self.k {
                    self.cur_pos = None;
                    None
                } else {
                    let buf_slice = &self.buffer[pos..pos + self.k];
                    
                    // advance for next time
                    if self.multiline && buf_slice.contains(&('\n' as u8)) {
                        // FIXME: this logic fails if kmer == 1; we repeat a character (or crash)
                        let mut buf: Vec<u8> = Vec::with_capacity(self.k);

                        for n in self.buffer[pos..].iter() {
                            if n != &('\n' as u8) && n != &('\r' as u8) {
                                buf.push(*n);
                            }
                            if buf.len() == self.k {
                                break;
                            }
                        }
                        
                        // FIXME: this may not work for a "\r\n" because it's only going ahead one?
                        // if the next position is going to start with whitespace, skip ahead
                        if self.buffer[pos + 1] == '\n' as u8 || self.buffer[pos + 1] == '\r' as u8 {
                            pos += 1
                        }

                        if pos > self.buffer.len() - self.k {
                            self.cur_pos = None;
                            None
                        } else {
                            self.cur_pos = Some(pos + 1);
                            Some(Cow::Owned(buf))
                        }
                    } else {
                        self.cur_pos = Some(pos + 1);
                        Some(Cow::Borrowed(buf_slice))
                    }
                }
            }
        }
    }
}
