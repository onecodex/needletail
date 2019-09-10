use std::borrow::Cow;

use memchr::{memchr, memchr2};

use crate::bitkmer::BitNuclKmer;
use crate::kmer::{complement, CanonicalKmers, Kmers};

pub fn normalize(seq: &[u8], allow_iupac: bool) -> Option<Vec<u8>> {
    //! Transform a FASTX sequence into it's "normalized" form.
    //!
    //! The normalized form is:
    //!  - only AGCTN and possibly - (for gaps)
    //!  - strip out any whitespace or line endings
    //!  - lowercase versions of these are uppercased
    //!  - U is converted to T (make everything a DNA sequence)
    //!  - some other punctuation is converted to gaps
    //!  - IUPAC bases may be converted to N's depending on the parameter passed in
    //!  - everything else is considered a N
    let mut buf: Vec<u8> = Vec::with_capacity(seq.len());
    let mut changed: bool = false;

    for n in seq.iter() {
        let (new_char, char_changed) = match (*n, allow_iupac) {
            c @ (b'A', _)
            | c @ (b'C', _)
            | c @ (b'G', _)
            | c @ (b'T', _)
            | c @ (b'N', _)
            | c @ (b'-', _) => (c.0, false),
            (b'a', _) => (b'A', true),
            (b'c', _) => (b'C', true),
            (b'g', _) => (b'G', true),
            // normalize uridine to thymine
            (b't', _) | (b'u', _) | (b'U', _) => (b'T', true),
            // normalize gaps
            (b'.', _) | (b'~', _) => (b'-', true),
            // logic for IUPAC bases (a little messy)
            c @ (b'B', true)
            | c @ (b'D', true)
            | c @ (b'H', true)
            | c @ (b'V', true)
            | c @ (b'R', true)
            | c @ (b'Y', true)
            | c @ (b'S', true)
            | c @ (b'W', true)
            | c @ (b'K', true)
            | c @ (b'M', true) => (c.0, false),
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

/// Mask tabs in header lines to `|`s
pub fn mask_header_tabs(id: &[u8]) -> Option<Vec<u8>> {
    memchr(b'\t', id).map(|_| {
        id.iter()
            .map(|x| if *x == b'\t' { b'|' } else { *x })
            .collect()
    })
}

/// Convert bad UTF8 characters into �s
pub fn mask_header_utf8(id: &[u8]) -> Option<Vec<u8>> {
    // this may potentially change the length of the id; we should probably
    // be doing something trickier like converting
    match String::from_utf8_lossy(id) {
        Cow::Owned(s) => Some(s.into_bytes()),
        Cow::Borrowed(_) => None,
    }
}

pub struct SequenceRecord<'a> {
    pub id: Cow<'a, [u8]>,
    pub seq: Cow<'a, [u8]>,
    pub qual: Option<Cow<'a, [u8]>>,
}

impl<'a> SequenceRecord<'a> {
    pub fn new(id: Cow<'a, [u8]>, seq: Cow<'a, [u8]>, qual: Option<Cow<'a, [u8]>>) -> Self {
        // there has to be a better way to do this?
        let cleaned_seq = match seq.strip_returns() {
            Cow::Owned(s) => Cow::Owned(s),
            Cow::Borrowed(_) => seq,
        };
        SequenceRecord {
            id,
            seq: cleaned_seq,
            qual,
        }
    }
}

impl<'a> From<&'a [u8]> for SequenceRecord<'a> {
    fn from(slice: &'a [u8]) -> Self {
        SequenceRecord::new(Cow::from(&b""[..]), slice.into(), None)
    }
}

/// A generic FASTX record that also abstracts over several logical operations
/// that can be performed on nucleic acid sequences.
pub trait Sequence<'a> {
    fn sequence(&'a self) -> &'a [u8];

    /// remove newlines from within FASTX records; currently the rate limiting step
    /// in FASTX parsing (in general; readfq also exhibits this behavior)
    fn strip_returns(&'a self) -> Cow<'a, [u8]> {
        let seq = self.sequence();

        // first part is a fast check to see if we need to do any allocations
        let mut i;
        match memchr2(b'\r', b'\n', &seq) {
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

    fn reverse_complement(&'a self) -> Vec<u8> {
        self.sequence()
            .iter()
            .rev()
            .map(|n| complement(*n))
            .collect()
    }

    fn normalize(&'a self, iupac: bool) -> Cow<'a, [u8]> {
        if let Some(s) = normalize(&self.sequence(), iupac) {
            s.into()
        } else {
            self.sequence().into()
        }
    }

    fn canonical_kmers(&'a self, k: u8, reverse_complement: &'a [u8]) -> CanonicalKmers<'a> {
        CanonicalKmers::new(self.sequence().as_ref(), reverse_complement, k)
    }

    fn kmers(&'a self, k: u8) -> Kmers<'a> {
        Kmers::new(self.sequence().as_ref(), k)
    }

    /// Return an iterator the returns valid kmers in 4-bit form
    fn bit_kmers(&'a self, k: u8, canonical: bool) -> BitNuclKmer<'a> {
        BitNuclKmer::new(self.sequence(), k, canonical)
    }
}

impl<'a> Sequence<'a> for &'a [u8] {
    fn sequence(&'a self) -> &'a [u8] {
        &self
    }
}

impl<'a> Sequence<'a> for [u8] {
    fn sequence(&'a self) -> &'a [u8] {
        &self
    }
}

impl<'a> Sequence<'a> for Cow<'a, [u8]> {
    fn sequence(&'a self) -> &'a [u8] {
        &self
    }
}

impl<'a> Sequence<'a> for SequenceRecord<'a> {
    fn sequence(&'a self) -> &'a [u8] {
        self.seq.as_ref()
    }
}

pub trait QualitySequence<'a>: Sequence<'a> {
    fn quality(&'a self) -> &'a [u8];

    /// Given a SeqRecord and a quality cutoff, mask out low-quality bases with
    /// `N` characters.
    ///
    /// Experimental.
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
        &self.0
    }
}

impl<'a> QualitySequence<'a> for (&'a [u8], &'a [u8]) {
    fn quality(&'a self) -> &'a [u8] {
        &self.1
    }
}

static EMPTY_VEC: &[u8] = b"";

impl<'a> QualitySequence<'a> for SequenceRecord<'a> {
    fn quality(&'a self) -> &'a [u8] {
        if let Some(q) = self.qual.as_ref() {
            q.as_ref()
        } else {
            &EMPTY_VEC
        }
        // fake high quality scores? vec![b'I'; self.sequence().len()]
    }
}

//
//    /// Fixes up potential problems with sequence headers including tabs being
//    /// present (may break downstream analyses with headers in TSVs) and with
//    /// non-UTF8 characters being present, e.g. non-breaking spaces on Windows
//    /// encodings (0x0A) breaks some tools.
//    pub fn mask_header(mut self) -> Self {
//        if let Some(id) = mask_header_tabs(&self.id) {
//            self.id = id.into();
//        }
//        if let Some(id) = mask_header_utf8(&self.id) {
//            self.id = id.into();
//        }
//        self
//    }

#[test]
fn test_quality_mask() {
    let seq_rec = (&b"AGCT"[..], &b"AAA0"[..]);
    let filtered_rec = seq_rec.quality_mask(b'5');
    assert_eq!(&filtered_rec[..], &b"AGCN"[..]);
}

#[test]
fn can_kmerize() {
    // test general function
    for (i, k) in b"AGCT".kmers(1).enumerate() {
        match i {
            0 => assert_eq!(k, &b"A"[..]),
            1 => assert_eq!(k, &b"G"[..]),
            2 => assert_eq!(k, &b"C"[..]),
            3 => assert_eq!(k, &b"T"[..]),
            _ => unreachable!("Too many kmers"),
        }
    }

    // test that we handle length 2 (and don't drop Ns)
    for (i, k) in b"ACNGT".kmers(2).enumerate() {
        match i {
            0 => assert_eq!(k, &b"AC"[..]),
            1 => assert_eq!(k, &b"CN"[..]),
            2 => assert_eq!(k, &b"NG"[..]),
            3 => assert_eq!(k, &b"GT"[..]),
            _ => unreachable!("Too many kmers"),
        }
    }

    // test that the minimum length works
    for k in b"AC".kmers(2) {
        assert_eq!(k, &b"AC"[..]);
    }
}

#[test]
fn can_canonicalize() {
    // test general function
    let seq = b"AGCT";
    for (i, (_, k, is_c)) in seq
        .canonical_kmers(1, &seq.reverse_complement())
        .enumerate()
    {
        match i {
            0 => {
                assert_eq!(k, &b"A"[..]);
                assert_eq!(is_c, false);
            }
            1 => {
                assert_eq!(k, &b"C"[..]);
                assert_eq!(is_c, true);
            }
            2 => {
                assert_eq!(k, &b"C"[..]);
                assert_eq!(is_c, false);
            }
            3 => {
                assert_eq!(k, &b"A"[..]);
                assert_eq!(is_c, true);
            }
            _ => unreachable!("Too many kmers"),
        }
    }

    let seq = b"AGCTA";
    for (i, (_, k, _)) in seq
        .canonical_kmers(2, &seq.reverse_complement())
        .enumerate()
    {
        match i {
            0 => assert_eq!(k, &b"AG"[..]),
            1 => assert_eq!(k, &b"GC"[..]),
            2 => assert_eq!(k, &b"AG"[..]),
            3 => assert_eq!(k, &b"TA"[..]),
            _ => unreachable!("Too many kmers"),
        }
    }

    let seq = b"AGNTA";
    for (i, (ix, k, _)) in seq
        .canonical_kmers(2, &seq.reverse_complement())
        .enumerate()
    {
        match i {
            0 => {
                assert_eq!(ix, 0);
                assert_eq!(k, &b"AG"[..]);
            }
            1 => {
                assert_eq!(ix, 3);
                assert_eq!(k, &b"TA"[..]);
            }
            _ => unreachable!("Too many kmers"),
        }
    }
}
