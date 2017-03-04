pub type BitKmerSeq = u64;
pub type BitKmer = (BitKmerSeq, u8);

/// Takes a BitKmer and adds a new base on the end, optionally loping off the
/// first base if the resulting kmer is too long.
fn extend_kmer(kmer: &mut BitKmer, new_char: &u8) -> bool {
    let new_char_int;
    match new_char {
        &b'A' | &b'a' => new_char_int = 0 as BitKmerSeq,
        &b'C' | &b'c' => new_char_int = 1 as BitKmerSeq,
        &b'G' | &b'g' => new_char_int = 2 as BitKmerSeq,
        &b'T' | &b't' => new_char_int = 3 as BitKmerSeq,
        _ => return false,
    };
    let new_kmer = (kmer.0 << 2) + new_char_int;

    // mask out any overflowed bits
    kmer.0 = new_kmer & (BitKmerSeq::pow(2, (2 * kmer.1) as u32) - 1) as BitKmerSeq;
    true
}

fn update_position(start_pos: &mut usize, kmer: &mut BitKmer, buffer: &[u8], initial: bool) -> bool {
    // check if we have enough "physical" space for one more kmer
    if *start_pos + kmer.1 as usize > buffer.len() {
        return false;
    }
    
    let mut kmer_len = (kmer.1 - 1) as usize;
    let mut stop_len = kmer.1 as usize;
    if initial {
        kmer_len = 0;
        stop_len = (kmer.1 - 1) as usize;
    }

    let mut cur_kmer = kmer;
    while kmer_len < stop_len {
        if extend_kmer(&mut cur_kmer, &buffer[*start_pos + kmer_len]) {
            kmer_len += 1;
        } else {
            kmer_len = 0;
            *cur_kmer = (0u64, cur_kmer.1);
            *start_pos += kmer_len + 1;
            if *start_pos + cur_kmer.1 as usize > buffer.len() {
                return false;
            }
        }
    }
    true
}


pub struct BitNuclKmer<'a> {
    start_pos: usize,
    cur_kmer: BitKmer,
    buffer: &'a [u8],
    canonical: bool,
}

impl<'a> BitNuclKmer<'a> {
    pub fn new(slice: &'a [u8], k: u8, canonical: bool) -> BitNuclKmer<'a> {
        let mut kmer = (0u64, k);
        let mut start_pos = 0;
        update_position(&mut start_pos, &mut kmer, slice, true);

        BitNuclKmer {
            start_pos: start_pos,
            cur_kmer: kmer,
            buffer: slice,
            canonical: canonical,
        }
    }
}

impl<'a> Iterator for BitNuclKmer<'a> {
    type Item = (usize, BitKmer,  bool);

    fn next(&mut self) -> Option<(usize, BitKmer, bool)> {
        if !update_position(&mut self.start_pos, &mut self.cur_kmer, self.buffer, false) {
            return None;
        }
        self.start_pos += 1;
        if self.canonical {
            let (kmer, was_rc) = canonical(self.cur_kmer);
            Some((self.start_pos - 1, kmer, was_rc))
        } else {
            Some((self.start_pos - 1, self.cur_kmer, false))
        }
    }
}

#[test]
fn can_kmerize() {
    // test general function
    let mut i = 0;
    for (_, k, _) in BitNuclKmer::new(b"AGCT", 1, false) {
        match i {
            0 => assert_eq!(k.0, 0b00 as BitKmerSeq),
            1 => assert_eq!(k.0, 0b10 as BitKmerSeq),
            2 => assert_eq!(k.0, 0b01 as BitKmerSeq),
            3 => assert_eq!(k.0, 0b11 as BitKmerSeq),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that we skip over N's
    i = 0;
    for (_, k, _) in BitNuclKmer::new(b"ACNGT", 2, false) {
        match i {
            0 => assert_eq!(k.0, 0b0001 as BitKmerSeq),
            1 => assert_eq!(k.0, 0b1011 as BitKmerSeq),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that we skip over N's and handle short kmers
    i = 0;
    for (_, k, _) in BitNuclKmer::new(b"ACNG", 2, false) {
        match i {
            0 => assert_eq!(k.0, 0x0001 as BitKmerSeq),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that the minimum length works
    i = 0;
    for (_, k, _) in BitNuclKmer::new(b"AC", 2, false) {
        match i {
            0 => assert_eq!(k.0, 0x0001 as BitKmerSeq),
            _ => assert!(false),
        }
        i += 1;
    }
}

#[test]
fn test_iterator() {
    let seq = "ACGTA".as_bytes();
    let mut kmer_iter = BitNuclKmer::new(seq, 3, false);
    assert_eq!(kmer_iter.next(), Some((0, (6, 3), false)));
    assert_eq!(kmer_iter.next(), Some((1, (27, 3), false)));
    assert_eq!(kmer_iter.next(), Some((2, (44, 3), false)));
    assert_eq!(kmer_iter.next(), None);

    let seq = "TA".as_bytes();
    let mut kmer_iter = BitNuclKmer::new(seq, 3, false);
    assert_eq!(kmer_iter.next(), None);
}


/// Reverse complement a BitKmer (reverses the sequence and swaps A<>T and G<>C)
pub fn reverse_complement(kmer: BitKmer) -> BitKmer {
    // FIXME: this is not going to work with BitKmers of u128 or u32
    // inspired from https://www.biostars.org/p/113640/
    let mut new_kmer = kmer.0;
    // reverse it
    new_kmer = (new_kmer >> 2 & 0x3333333333333333) | (new_kmer & 0x3333333333333333) << 2;
    new_kmer = (new_kmer >> 4 & 0x0F0F0F0F0F0F0F0F) | (new_kmer & 0x0F0F0F0F0F0F0F0F) << 4;
    new_kmer = (new_kmer >> 8 & 0x00FF00FF00FF00FF) | (new_kmer & 0x00FF00FF00FF00FF) << 8;
    new_kmer = (new_kmer >> 16 & 0x0000FFFF0000FFFF) | (new_kmer & 0x0000FFFF0000FFFF) << 16;
    new_kmer = (new_kmer >> 32 & 0x00000000FFFFFFFF) | (new_kmer & 0x00000000FFFFFFFF) << 32;
    // complement it
    new_kmer ^= 0xFFFFFFFFFFFFFFFF;
    // shift it to the right size
    new_kmer = new_kmer >> (2 * (32 - kmer.1));
    (new_kmer, kmer.1)
}

#[test]
fn test_reverse_complement() {
  assert_eq!(reverse_complement((0b000000, 3)).0, 0b111111);
  assert_eq!(reverse_complement((0b111111, 3)).0, 0b000000);
  assert_eq!(reverse_complement((0b00000000, 4)).0, 0b11111111);
  assert_eq!(reverse_complement((0b00011011, 4)).0, 0b00011011);
}

/// Return the lexigraphically lowest of the BitKmer and its reverse complement and 
/// whether the returned kmer is the reverse_complement (true) or the original (false)
pub fn canonical(kmer: BitKmer) -> (BitKmer, bool) {
    let rc = reverse_complement(kmer);
    if kmer.0 > rc.0 {
        (rc, true)
    } else {
        (kmer, false)
    }
}

/// Find the lexigraphically lowest substring of a given length in the BitKmer
pub fn minimizer(kmer: BitKmer, minmer_size: u8) -> BitKmer {
    let mut new_kmer = kmer.0;
    let mut lowest = !(0 as BitKmerSeq);
    let bitmask = (BitKmerSeq::pow(2, (2 * minmer_size) as u32) - 1) as BitKmerSeq;
    for _ in 0..(kmer.1 - minmer_size + 1) {
        let cur = bitmask & new_kmer;
        if cur < lowest {
            lowest = cur;
        }
        let cur_rev = reverse_complement((bitmask & new_kmer, kmer.1));
        if cur_rev.0 < lowest {
            lowest = cur_rev.0;
        }
        new_kmer >>= 2;
    }
    (lowest, kmer.1)
}

#[test]
fn test_minimizer() {
    assert_eq!(minimizer((0b001011, 3), 2).0, 0b0010);
    assert_eq!(minimizer((0b001011, 3), 1).0, 0b00);
    assert_eq!(minimizer((0b11000011, 4), 2).0, 0b0000);
    assert_eq!(minimizer((0b110001, 3), 2).0, 0b0001);
}

pub fn bitmer_to_bytes(kmer: BitKmer) -> Vec<u8> {
    let mut new_kmer = kmer.0;
    let mut new_kmer_str = Vec::new();
    // we're reading the bases off from the "high" end of the integer so we need to do some
    // math to figure out where they start (this helps us just pop the bases on the end
    // of the working buffer as we read them off "left to right")
    let offset = (kmer.1 - 1) * 2;
    let bitmask = BitKmerSeq::pow(2, (2 * kmer.1 - 1) as u32) + BitKmerSeq::pow(2, (2 * kmer.1 - 2) as u32);

    for _ in 0..kmer.1 {
        let new_char = (new_kmer & bitmask) >> offset;
        new_kmer <<= 2;
        new_kmer_str.push(match new_char {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => panic!("Mathematical impossibility"),
        });
    }
    new_kmer_str
}

#[test]
fn test_bitmer_to_bytes() {
    assert_eq!(bitmer_to_bytes((1 as BitKmerSeq, 1)), b"C".to_vec());
    assert_eq!(bitmer_to_bytes((60 as BitKmerSeq, 3)), b"TTA".to_vec());
    assert_eq!(bitmer_to_bytes((0 as BitKmerSeq, 3)), b"AAA".to_vec());
}

pub fn bytes_to_bitmer(kmer: &[u8]) -> BitKmer {
    let k = kmer.len() as u8;

    let mut bit_kmer = (0u64, k);
    for i in 0..k {
        extend_kmer(&mut bit_kmer, &kmer[i as usize]);
    }
    bit_kmer
}

#[test]
fn test_bytes_to_bitkmer() {
    let mut ikmer: BitKmer = bytes_to_bitmer("C".as_bytes());
    assert_eq!(ikmer.0, 1 as BitKmerSeq);

    ikmer = bytes_to_bitmer("TTA".as_bytes());
    assert_eq!(ikmer.0, 60 as BitKmerSeq);

    ikmer = bytes_to_bitmer("AAA".as_bytes());
    assert_eq!(ikmer.0, 0 as BitKmerSeq);
}
