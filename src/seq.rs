use std::borrow::Cow;

use bitkmer::BitNuclKmer;
use kmer::{is_good_base, complement, normalize};

/// A generic FASTX record that also abstracts over several logical operations
/// that can be performed on nucleic acid sequences.
pub struct SeqRecord<'a> {
    pub id: &'a str,
    pub seq: Cow<'a, [u8]>,
    pub qual: Option<&'a [u8]>,
    rev_seq: Option<Vec<u8>>,
}

impl<'a> SeqRecord<'a> {
    pub fn new(id: &'a str, seq: Cow<'a, [u8]>, qual: Option<&'a [u8]>) -> Self {
        SeqRecord { id: id, seq: seq, qual: qual, rev_seq: None }
    }

    pub fn from_bytes(seq: &'a [u8]) -> Self {
        SeqRecord { id: "", seq: Cow::Borrowed(seq), qual: None, rev_seq: None }
    }

    /// Given a SeqRecord and a quality cutoff, mask out low-quality bases with
    /// `N` characters.
    ///
    /// Experimental.
    pub fn quality_mask(self, ref score: u8) -> Self {
        match self.qual {
            None => self,
            Some(quality) => {
                // could maybe speed this up by doing a copy of base and then
                // iterating though qual and masking?
                let seq = self.seq
                    .iter()
                    .zip(quality.iter())
                    .map(|(base, qual)| {
                        if qual < score {
                            b'N'
                        } else {
                            base.clone()
                        }
                    })
                    .collect();
                SeqRecord {
                    id: self.id,
                    seq: seq,
                    qual: self.qual,
                    rev_seq: None,
                }
            },
        }
    }

    /// Capitalize everything and mask unknown bases to N
    pub fn normalize(self, iupac: bool) -> Self {
        let seq = normalize(&self.seq, iupac);
        SeqRecord {
            id: self.id,
            seq: Cow::Owned(seq),
            qual: self.qual,
            rev_seq: None,
        }
    }

    /// Return an iterator the returns valid kmers
    pub fn kmers<'b, 'c>(&'b mut self, k: u8, canonical: bool) -> NuclKmer<'c> where 'b: 'c {
        if canonical {
            self.rev_seq = Some(self.seq.iter().rev().map(|n| complement(n)).collect());
        }
        match self.rev_seq {
            Some(ref rev_seq) => NuclKmer::new(&self.seq, Some(&rev_seq), k),
            None => NuclKmer::new(&self.seq, None, k),
        }
    }

    /// Return an iterator the returns valid kmers in 4-bit form
    pub fn bit_kmers<'b>(&'b self, k: u8, canonical: bool) -> BitNuclKmer<'b> {
        BitNuclKmer::new(&self.seq, k, canonical)
    }
}


pub struct NuclKmer<'a> {
    k: u8,
    start_pos: usize,
    buffer: &'a [u8],
    rc_buffer: Option<&'a [u8]>,
}

fn update_position(start_pos: &mut usize, k: u8, buffer: &[u8], initial: bool) -> bool {
    // check if we have enough "physical" space for one more kmer
    if *start_pos + k as usize > buffer.len() {
        return false;
    }

    let mut kmer_len = (k - 1) as usize;
    let mut stop_len = k as usize;
    if initial {
        kmer_len = 0;
        stop_len = (k - 1) as usize;
    }

    while kmer_len < stop_len {
        if is_good_base(buffer[*start_pos + kmer_len]) {
            kmer_len += 1;
        } else {
            kmer_len = 0;
            *start_pos += kmer_len + 1;
            if *start_pos + k as usize > buffer.len() {
                return false;
            }
        }
    }
    true
}

impl<'a> NuclKmer<'a> {
    //! A kmer-izer for a nucleotide/amino acid sequence; returning slices to the original data
    pub fn new(buffer: &'a [u8], rc_buffer: Option<&'a [u8]>, k: u8) -> NuclKmer<'a> {
        let mut start_pos = 0;
        update_position(&mut start_pos, k, buffer, true);
        NuclKmer {
            k: k,
            start_pos: start_pos,
            buffer: buffer,
            rc_buffer: rc_buffer,
        }
    }
}

impl<'a> Iterator for NuclKmer<'a> {
    type Item = (usize, &'a [u8], bool);

    fn next(&mut self) -> Option<(usize, &'a [u8], bool)> {
        if !update_position(&mut self.start_pos, self.k, self.buffer, false) {
            return None;
        }
        let pos = self.start_pos;
        self.start_pos += 1;

        let result = &self.buffer[pos..pos + self.k as usize];
        match self.rc_buffer {
            None => Some((pos, result, false)),
            Some(rc_buffer) => {
                let rc_result = &rc_buffer[rc_buffer.len() - pos - self.k as usize..rc_buffer.len() - pos];
                if result < rc_result {
                    Some((pos, result, false))
                } else {
                    Some((pos, rc_result, true))
                }
            }
        }
    }
}

#[test]
fn test_quality_mask() {
    let seq_rec = SeqRecord {
        id: "",
        seq: Cow::Borrowed(&b"AGCT"[..]),
        qual: Some(&b"AAA0"[..]),
        rev_seq: None,
    };
    let filtered_rec = seq_rec.quality_mask('5' as u8);
    assert_eq!(&filtered_rec.seq[..], &b"AGCN"[..]);
}

#[test]
fn can_kmerize() {
    // test general function
    let mut i = 0;
    for (_, k, _) in SeqRecord::from_bytes(b"AGCT").kmers(1, false) {
        match i {
            0 => assert_eq!(k, &b"A"[..]),
            1 => assert_eq!(k, &b"G"[..]),
            2 => assert_eq!(k, &b"C"[..]),
            3 => assert_eq!(k, &b"T"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that we skip over N's
    i = 0;
    for (_, k, _) in SeqRecord::from_bytes(b"ACNGT").kmers(2, false) {
        match i {
            0 => assert_eq!(k, &b"AC"[..]),
            1 => assert_eq!(k, &b"GT"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    // test that we skip over N's and handle short kmers
    i = 0;
    for (ix, k, _) in SeqRecord::from_bytes(b"ACNG").kmers(2, false) {
        match i {
            0 => {
                assert_eq!(ix, 0);
                assert_eq!(k, &b"AC"[..]);
            },
            _ => assert!(false),
        }
        i += 1;
    }

    // test that the minimum length works
    for (_, k, _) in SeqRecord::from_bytes(b"AC").kmers(2, false) {
        assert_eq!(k, &b"AC"[..]);
    }
}

#[test]
fn can_canonicalize() {
    // test general function
    let mut i = 0;
    for (_, k, is_c) in SeqRecord::from_bytes(b"AGCT").kmers(1, true) {
        match i {
            0 => {
                assert_eq!(k, &b"A"[..]);
                assert_eq!(is_c, false);
            },
            1 => {
                assert_eq!(k, &b"C"[..]);
                assert_eq!(is_c, true);
            },
            2 => {
                assert_eq!(k, &b"C"[..]);
                assert_eq!(is_c, false);
            },
            3 => {
                assert_eq!(k, &b"A"[..]);
                assert_eq!(is_c, true);
            },
            _ => assert!(false),
        }
        i += 1;
    }

    let mut i = 0;
    for (_, k, _) in SeqRecord::from_bytes(b"AGCTA").kmers(2, true) {
        match i {
            0 => assert_eq!(k, &b"AG"[..]),
            1 => assert_eq!(k, &b"GC"[..]),
            2 => assert_eq!(k, &b"AG"[..]),
            3 => assert_eq!(k, &b"TA"[..]),
            _ => assert!(false),
        }
        i += 1;
    }

    let mut i = 0;
    for (ix, k, _) in SeqRecord::from_bytes(b"AGNTA").kmers(2, true) {
        match i {
            0 => {
                assert_eq!(ix, 0);
                assert_eq!(k, &b"AG"[..]);
            },
            1 => {
                assert_eq!(ix, 3);
                assert_eq!(k, &b"TA"[..]);
            },
            _ => assert!(false),
        }
        i += 1;
    }
}
