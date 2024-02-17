//! Generic functions for working with (primarily nucleic acid) sequences
use std::borrow::Cow;

use memchr::memchr2;

use crate::bitkmer::BitNuclKmer;
use crate::kmer::{CanonicalKmers, Kmers};

/// Transform a nucleic acid sequence into its "normalized" form.
///
/// The normalized form is:
///  - only AGCTN and possibly - (for gaps)
///  - strip out any whitespace or line endings
///  - lowercase versions of these are uppercased
///  - U is converted to T (make everything a DNA sequence)
///  - some other punctuation is converted to gaps
///  - IUPAC bases may be converted to N's depending on the parameter passed in
///  - everything else is considered a N
pub fn normalize(seq: &[u8], allow_iupac: bool) -> Option<Vec<u8>> {
    let mut buf: Vec<u8> = Vec::with_capacity(seq.len());
    let mut changed: bool = false;

    for n in seq.iter() {
        let (new_char, char_changed) = match (*n, allow_iupac) {
            c @ (b'A' | b'C' | b'G' | b'T' | b'N' | b'-', _) => (c.0, false),
            (b'a', _) => (b'A', true),
            (b'c', _) => (b'C', true),
            (b'g', _) => (b'G', true),
            // normalize uridine to thymine
            (b't' | b'u' | b'U', _) => (b'T', true),
            // normalize gaps
            (b'.' | b'~', _) => (b'-', true),
            // logic for IUPAC bases (a little messy)
            c @ (b'B' | b'D' | b'H' | b'V' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M', true) => {
                (c.0, false)
            }
            (b'b', true) => (b'B', true),
            (b'd', true) => (b'D', true),
            (b'h', true) => (b'H', true),
            (b'v', true) => (b'V', true),
            (b'r', true) => (b'R', true),
            (b'y', true) => (b'Y', true),
            (b's', true) => (b'S', true),
            (b'w', true) => (b'W', true),
            (b'k', true) => (b'K', true),
            (b'm', true) => (b'M', true),
            // remove all whitespace and line endings
            (b' ', _) | (b'\t', _) | (b'\r', _) | (b'\n', _) => (b' ', true),
            // everything else is an N
            _ => (b'N', true),
        };
        changed = changed || char_changed;
        if new_char != b' ' {
            buf.push(new_char);
        }
    }
    if changed {
        Some(buf)
    } else {
        None
    }
}

/// Returns the complementary base for a given IUPAC base code.
///
/// Does not work for RNA sequences (maybe we should raise an error or something?)
#[inline]
pub fn complement(n: u8) -> u8 {
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

/// Taking in a sequence string, return the canonical form of the sequence
/// (e.g. the lexigraphically lowest of either the original sequence or its
/// reverse complement)
pub fn canonical(seq: &[u8]) -> Cow<[u8]> {
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

/// A generic FASTX record that also abstracts over several logical operations
/// that can be performed on nucleic acid sequences.
pub trait Sequence<'a> {
    fn sequence(&'a self) -> &'a [u8];

    /// Remove newlines from the sequence; this handles `\r`, `\n`, and `\r\n`
    /// and removes internal newlines in addition to ones at the end.
    /// Primarily used for FASTA multiline records, but can also help process
    /// (the much rarer) multiline FASTQs. Always use before iteration methods
    /// below to ensure no newlines are being returned with e.g. `.kmers`.
    /// If you are using `normalize`, you do not need to call this function directly.
    fn strip_returns(&'a self) -> Cow<'a, [u8]> {
        let seq = self.sequence();

        // first part is a fast check to see if we need to do any allocations
        let mut i;
        match memchr2(b'\r', b'\n', seq) {
            Some(break_loc) => i = break_loc,
            None => return seq.into(),
        }
        // we found a newline; create a new buffer and stripping out newlines
        // and writing into it
        let mut new_buf = Vec::with_capacity(seq.len() - 1);
        new_buf.extend_from_slice(&seq[..i]);
        while i < seq.len() {
            match memchr2(b'\r', b'\n', &seq[i..]) {
                None => {
                    new_buf.extend_from_slice(&seq[i..]);
                    break;
                }
                Some(match_pos) => {
                    new_buf.extend_from_slice(&seq[i..i + match_pos]);
                    i += match_pos + 1;
                }
            }
        }
        new_buf.into()
    }

    /// Returns the reverse complement of a sequence. Biologically this is
    /// equivalent to the sequence of the strand opposite the one you pass
    /// in.
    ///
    /// ```
    /// use needletail::Sequence;
    ///
    /// assert_eq!(b"AACC".reverse_complement(), b"GGTT");
    /// ```
    fn reverse_complement(&'a self) -> Vec<u8> {
        self.sequence()
            .iter()
            .rev()
            .map(|n| complement(*n))
            .collect()
    }

    /// [Nucleic Acids] Normalizes the sequence. See documentation for
    /// `needletail::sequence::normalize`. Do not use on amino acid
    /// sequences. Note that this returns a Cow so you may have to coerce
    /// to a Vec<u8> or &[u8] as necessary.
    ///
    /// ```
    /// use needletail::Sequence;
    ///
    /// // IUPAC bases are coerced to N's if `false`
    /// assert_eq!(b"ADGH".normalize(false).as_ref(), b"ANGN");
    /// // otherwise they're preserved
    /// assert_eq!(b"ADGH".normalize(true).as_ref(), b"ADGH");
    ///
    /// // Uridine residues are converted to thymidine
    /// assert_eq!(b"ACGU".normalize(true).as_ref(), b"ACGT");
    /// ```
    fn normalize(&'a self, iupac: bool) -> Cow<'a, [u8]> {
        if let Some(s) = normalize(self.sequence(), iupac) {
            s.into()
        } else {
            self.sequence().into()
        }
    }

    /// [Nucleic Acids] Returns an iterator over the sequence that skips
    /// non-ACGT bases and returns a tuple containing (position, the
    /// canonicalized kmer, if the sequence is the complement of the original).
    fn canonical_kmers(&'a self, k: u8, reverse_complement: &'a [u8]) -> CanonicalKmers<'a> {
        CanonicalKmers::new(self.sequence(), reverse_complement, k)
    }

    /// Returns an iterator that returns a sliding window of k-sized
    /// sequences (k-mers). Does not skip whitespace or correct bases in the
    /// original sequence so `.normalize` or `.strip_returns` may be
    /// appropriate to use first.
    fn kmers(&'a self, k: u8) -> Kmers<'a> {
        Kmers::new(self.sequence(), k)
    }

    /// Return an iterator that returns valid kmers in 4-bit form
    fn bit_kmers(&'a self, k: u8, canonical: bool) -> BitNuclKmer<'a> {
        BitNuclKmer::new(self.sequence(), k, canonical)
    }
}

impl<'a> Sequence<'a> for &'a [u8] {
    fn sequence(&'a self) -> &'a [u8] {
        self
    }
}

impl<'a> Sequence<'a> for [u8] {
    fn sequence(&'a self) -> &'a [u8] {
        self
    }
}

impl<'a> Sequence<'a> for Cow<'a, [u8]> {
    fn sequence(&'a self) -> &'a [u8] {
        self
    }
}

/// [⚠️Unstable] A trait to wrap over sequence data that has associated
/// quality information.
///
/// Will be stabilized once we figure out a good way to handle sequences that
/// have _optional_ quality information (like SequenceRecord) because the
/// return trait requires a slice from an immutable reference and
/// SequenceRecords can't return that without modifying themselves.
pub trait QualitySequence<'a>: Sequence<'a> {
    fn quality(&'a self) -> &'a [u8];

    /// Given a SeqRecord and a quality cutoff, mask out low-quality bases with
    /// `N` characters.
    fn quality_mask(&'a self, score: u8) -> Cow<'a, [u8]> {
        let qual = self.quality();
        // could maybe speed this up by doing a copy of base and then
        // iterating though qual and masking?
        let seq: Vec<u8> = self
            .sequence()
            .iter()
            .zip(qual.iter())
            .map(|(base, qual)| if *qual < score { b'N' } else { *base })
            .collect();
        seq.into()
    }
}

impl<'a> Sequence<'a> for (&'a [u8], &'a [u8]) {
    fn sequence(&'a self) -> &'a [u8] {
        self.0
    }
}

impl<'a> QualitySequence<'a> for (&'a [u8], &'a [u8]) {
    fn quality(&'a self) -> &'a [u8] {
        self.1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_normalize() {
        assert_eq!(normalize(b"ACGTU", false), Some(b"ACGTT".to_vec()));
        assert_eq!(normalize(b"acgtu", false), Some(b"ACGTT".to_vec()));

        assert_eq!(normalize(b"N.N-N~N N", false), Some(b"N-N-N-NN".to_vec()));

        assert_eq!(normalize(b"BDHVRYSWKM", true), None);
        assert_eq!(normalize(b"bdhvryswkm", true), Some(b"BDHVRYSWKM".to_vec()));
        assert_eq!(
            normalize(b"BDHVRYSWKM", false),
            Some(b"NNNNNNNNNN".to_vec())
        );
        assert_eq!(
            normalize(b"bdhvryswkm", false),
            Some(b"NNNNNNNNNN".to_vec())
        );
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'a'), b't');
        assert_eq!(complement(b'c'), b'g');
        assert_eq!(complement(b'g'), b'c');
        assert_eq!(complement(b'n'), b'n');
    }

    #[test]
    fn can_canonicalize() {
        assert_eq!(canonical(b"A"), Cow::Borrowed(b"A"));
        assert_eq!(canonical(b"T"), Cow::Owned::<[u8]>(b"A".to_vec()));
        assert_eq!(canonical(b"AAGT"), Cow::Borrowed(b"AAGT"));
        assert_eq!(canonical(b"ACTT"), Cow::Owned::<[u8]>(b"AAGT".to_vec()));
        assert_eq!(canonical(b"GC"), Cow::Borrowed(b"GC"));
    }

    #[test]
    fn can_minimize() {
        let minmer = minimizer(&b"ATTTCG"[..], 3);
        assert_eq!(&minmer[..], b"AAA");
    }

    #[test]
    fn test_quality_mask() {
        let seq_rec = (&b"AGCT"[..], &b"AAA0"[..]);
        let filtered_rec = seq_rec.quality_mask(b'5');
        assert_eq!(&filtered_rec[..], &b"AGCN"[..]);
    }
}
