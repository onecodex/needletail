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

#[test]
fn test_complement() {
    assert_eq!(complement(b'a'), b't');
    assert_eq!(complement(b'c'), b'g');
    assert_eq!(complement(b'g'), b'c');
    assert_eq!(complement(b'n'), b'n');
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

#[test]
fn can_canonicalize() {
    assert!(canonical(b"A") == Cow::Borrowed(b"A"));
    assert!(canonical(b"T") == Cow::Owned::<[u8]>(b"A".to_vec()));
    assert!(canonical(b"AAGT") == Cow::Borrowed(b"AAGT"));
    assert!(canonical(b"ACTT") == Cow::Owned::<[u8]>(b"AAGT".to_vec()));
    assert!(canonical(b"GC") == Cow::Borrowed(b"GC"));
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

#[test]
fn can_minimize() {
    let minmer = minimizer(&b"ATTTCG"[..], 3);
    assert_eq!(&minmer[..], b"AAA");
}

pub fn is_good_base(chr: u8) -> bool {
    match chr as char {
        'a' | 'c' | 'g' | 't' | 'A' | 'C' | 'G' | 'T' => true,
        _ => false,
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

    let (mut kmer_len, stop_len) = if initial {
        (0, (k - 1) as usize)
    } else {
        ((k - 1) as usize, k as usize)
    };

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
            k,
            start_pos,
            buffer,
            rc_buffer,
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
                let rc_result =
                    &rc_buffer[rc_buffer.len() - pos - self.k as usize..rc_buffer.len() - pos];
                if result < rc_result {
                    Some((pos, result, false))
                } else {
                    Some((pos, rc_result, true))
                }
            },
        }
    }
}
