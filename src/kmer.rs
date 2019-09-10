//! This module contains functions for kmerizing a longer sequence and various
//! utilities for dealing with these kmers.
//!
//!
use std::borrow::Cow;

#[inline]
pub fn complement(n: u8) -> u8 {
    //! Returns the complementary base for a given IUPAC base code.
    //!
    //! Does not work for RNA sequences (maybe we should raise an error or something?)
    match n {
        b'a' => b't',
        b'A' => b'T',
        b'c' => b'g',
        b'C' => b'G',
        b'g' => b'c',
        b'G' => b'C',
        b't' => b'a',
        b'T' => b'A',

        // IUPAC codes
        b'r' => b'y',
        b'y' => b'r',
        b'k' => b'm',
        b'm' => b'k',
        b'b' => b'v',
        b'v' => b'b',
        b'd' => b'h',
        b'h' => b'd',
        b's' => b's',
        b'w' => b'w',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        b'S' => b'S',
        b'W' => b'W',

        // anything else just pass through
        // 'u' | 'U' => panic!("Does not support complements of U"),
        x => x,
    }
}

pub fn canonical(seq: &[u8]) -> Cow<[u8]> {
    //! Taking in a sequence string, return the canonical form of the sequence
    //! (e.g. the lexigraphically lowest of either the original sequence or its
    //! reverse complement)
    let mut buf: Vec<u8> = Vec::with_capacity(seq.len());
    // enough just keeps our comparisons from happening after they need to
    let mut enough = false;
    let mut original_was_canonical = false;

    // loop through the kmer and its reverse complement simultaneously
    for (rn, n) in seq.iter().rev().map(|n| complement(*n)).zip(seq.iter()) {
        buf.push(rn);
        if !enough && (*n < rn) {
            original_was_canonical = true;
            break;
        } else if !enough && (rn < *n) {
            enough = true;
        }
        // unstated if branch: if rn == n, keep comparing
    }
    match (original_was_canonical, enough) {
        (true, true) => panic!("Bug: should never set original_was_canonical if enough == true"),
        (true, false) => seq.into(),
        (false, true) => buf.into(),
        // the sequences were completely equal, return the ref
        (false, false) => seq.into(),
    }
}

/// Find the lexigraphically smallest substring of `seq` of length `length`
///
/// There's probably a faster algorithm for this somewhere...
pub fn minimizer(seq: &[u8], length: usize) -> Cow<[u8]> {
    let reverse_complement: Vec<u8> = seq.iter().rev().map(|n| complement(*n)).collect();
    let mut minmer = Cow::Borrowed(&seq[..length]);

    for (kmer, rc_kmer) in seq.windows(length).zip(reverse_complement.windows(length)) {
        if *kmer < minmer[..] {
            minmer = kmer.into();
        }
        if *rc_kmer < minmer[..] {
            minmer = rc_kmer.to_vec().into();
        }
    }
    minmer
}

pub fn is_good_base(chr: u8) -> bool {
    match chr as char {
        'a' | 'c' | 'g' | 't' | 'A' | 'C' | 'G' | 'T' => true,
        _ => false,
    }
}

pub struct Kmers<'a> {
    k: u8,
    start_pos: usize,
    buffer: &'a [u8],
}

impl<'a> Kmers<'a> {
    //! A kmer-izer for a nucleotide/amino acid sequence; returning slices to the original data
    pub fn new(buffer: &'a [u8], k: u8) -> Self {
        Kmers {
            k,
            start_pos: 0,
            buffer,
        }
    }
}

impl<'a> Iterator for Kmers<'a> {
    type Item = &'a [u8];

    fn next(&mut self) -> Option<Self::Item> {
        if self.start_pos + self.k as usize > self.buffer.len() {
            return None;
        }
        let pos = self.start_pos;
        self.start_pos += 1;
        Some(&self.buffer[pos..pos + self.k as usize])
    }
}

pub struct CanonicalKmers<'a> {
    k: u8,
    start_pos: usize,
    buffer: &'a [u8],
    rc_buffer: &'a [u8],
}

/// A kmer-izer for a nucleotide acid sequences to return canonical kmers.
/// Returns the position of the kmer, a slice to the original data, and
/// an boolean indicating if the kmer returned is the original or the reverse
/// complement.
impl<'a> CanonicalKmers<'a> {
    pub fn new(buffer: &'a [u8], rc_buffer: &'a [u8], k: u8) -> Self {
        let mut nucl_kmers = CanonicalKmers {
            k,
            start_pos: 0,
            buffer,
            rc_buffer,
        };
        nucl_kmers.update_position(true);
        nucl_kmers
    }

    fn update_position(&mut self, initial: bool) -> bool {
        // check if we have enough "physical" space for one more kmer
        if self.start_pos + self.k as usize > self.buffer.len() {
            return false;
        }

        let (mut kmer_len, stop_len) = if initial {
            (0, (self.k - 1) as usize)
        } else {
            ((self.k - 1) as usize, self.k as usize)
        };

        while kmer_len < stop_len {
            if is_good_base(self.buffer[self.start_pos + kmer_len]) {
                kmer_len += 1;
            } else {
                kmer_len = 0;
                self.start_pos += kmer_len + 1;
                if self.start_pos + self.k as usize > self.buffer.len() {
                    return false;
                }
            }
        }
        true
    }
}

impl<'a> Iterator for CanonicalKmers<'a> {
    type Item = (usize, &'a [u8], bool);

    fn next(&mut self) -> Option<(usize, &'a [u8], bool)> {
        if !self.update_position(false) {
            return None;
        }
        let pos = self.start_pos;
        self.start_pos += 1;

        let result = &self.buffer[pos..pos + self.k as usize];
        let rc_buffer = self.rc_buffer;
        let rc_result = &rc_buffer[rc_buffer.len() - pos - self.k as usize..rc_buffer.len() - pos];
        if result < rc_result {
            Some((pos, result, false))
        } else {
            Some((pos, rc_result, true))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'a'), b't');
        assert_eq!(complement(b'c'), b'g');
        assert_eq!(complement(b'g'), b'c');
        assert_eq!(complement(b'n'), b'n');
    }

    #[test]
    fn can_canonicalize() {
        assert!(canonical(b"A") == Cow::Borrowed(b"A"));
        assert!(canonical(b"T") == Cow::Owned::<[u8]>(b"A".to_vec()));
        assert!(canonical(b"AAGT") == Cow::Borrowed(b"AAGT"));
        assert!(canonical(b"ACTT") == Cow::Owned::<[u8]>(b"AAGT".to_vec()));
        assert!(canonical(b"GC") == Cow::Borrowed(b"GC"));
    }

    #[test]
    fn can_minimize() {
        let minmer = minimizer(&b"ATTTCG"[..], 3);
        assert_eq!(&minmer[..], b"AAA");
    }

}
