use std::borrow::Cow;

use memchr::memchr;

use crate::bitkmer::BitNuclKmer;
use crate::kmer::{complement, NuclKmer};

pub fn normalize(seq: &[u8], iupac: bool) -> (Vec<u8>, bool) {
    //! Transform a FASTX sequence into it's "normalized" form.
    //!
    //! The normalized form is:
    //!  - only AGCTN and possibly . (for gaps)
    //!  - lowercase versions of these are uppercased
    //!  - U is converted to T (make everything a DNA sequence)
    //!  - some other punctuation is converted to gaps
    //!  - IUPAC bases may be converted to N's depending on the parameter passed in
    //!  - everything else is considered a N
    let mut buf: Vec<u8> = Vec::with_capacity(seq.len());
    let mut changed: bool = false;

    for n in seq.iter() {
        let (new_char, char_changed) = match (*n, iupac) {
            c @ (b'A', _)
            | c @ (b'C', _)
            | c @ (b'G', _)
            | c @ (b'T', _)
            | c @ (b'N', _)
            | c @ (b'.', _) => (c.0, false),
            (b'a', _) => (b'A', true),
            (b'c', _) => (b'C', true),
            (b'g', _) => (b'G', true),
            // normalize uridine to thymine
            (b't', _) | (b'u', _) | (b'U', _) => (b'T', true),
            // normalize gaps
            (b'-', _) | (b'~', _) | (b' ', _) => (b'.', true),
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
            _ => (b'N', true),
        };
        changed = changed || char_changed;
        buf.push(new_char);
    }
    (buf, changed)
}

#[test]
fn test_normalize() {
    assert_eq!(normalize(b"ACGTU", false), (b"ACGTT".to_vec(), true));
    assert_eq!(normalize(b"acgtu", false), (b"ACGTT".to_vec(), true));

    assert_eq!(
        normalize(b"N.N-N~N N", false),
        (b"N.N.N.N.N".to_vec(), true)
    );

    assert_eq!(
        normalize(b"BDHVRYSWKM", true),
        (b"BDHVRYSWKM".to_vec(), false)
    );
    assert_eq!(
        normalize(b"bdhvryswkm", true),
        (b"BDHVRYSWKM".to_vec(), true)
    );
    assert_eq!(
        normalize(b"BDHVRYSWKM", false),
        (b"NNNNNNNNNN".to_vec(), true)
    );
    assert_eq!(
        normalize(b"bdhvryswkm", false),
        (b"NNNNNNNNNN".to_vec(), true)
    );
}

/// A generic FASTX record that also abstracts over several logical operations
/// that can be performed on nucleic acid sequences.
#[derive(Clone, Debug)]
pub struct SeqRecord<'a> {
    pub id: Cow<'a, str>,
    pub seq: Cow<'a, [u8]>,
    pub qual: Option<Cow<'a, [u8]>>,
    rev_seq: Option<Vec<u8>>,
}

impl<'a> SeqRecord<'a> {
    pub fn new(id: &'a str, seq: Cow<'a, [u8]>, qual: Option<&'a [u8]>) -> Self {
        SeqRecord {
            id: id.into(),
            seq,
            qual: qual.map(Cow::Borrowed),
            rev_seq: None,
        }
    }

    pub fn from_bytes(seq: &'a [u8]) -> Self {
        SeqRecord {
            id: "".into(),
            seq: seq.into(),
            qual: None,
            rev_seq: None,
        }
    }

    /// Given a SeqRecord and a quality cutoff, mask out low-quality bases with
    /// `N` characters.
    ///
    /// Experimental.
    pub fn quality_mask(self, score: u8) -> Self {
        if self.qual == None {
            return self;
        }
        let qual = self.qual.unwrap().into_owned();
        // could maybe speed this up by doing a copy of base and then
        // iterating though qual and masking?
        let seq = self
            .seq
            .iter()
            .zip(qual.iter())
            .map(|(base, qual)| if *qual < score { b'N' } else { *base })
            .collect();
        SeqRecord {
            id: self.id,
            seq,
            qual: Some(Cow::Owned(qual)),
            rev_seq: None,
        }
    }

    /// Capitalize everything and mask unknown bases to N
    pub fn normalize(mut self, iupac: bool) -> Self {
        let (seq, changed) = normalize(&self.seq, iupac);
        if changed {
            self.seq = seq.into();
        }
        self
    }

    /// Mask tabs in header lines to `|`s
    ///
    /// Returns `true` if the header was masked
    pub fn mask_header(mut self) -> Self {
        if memchr(b'\t', self.id.as_ref().as_bytes()).is_some() {
            self.id = self.id.as_ref().replace("\t", "|").into();
        }
        self
    }

    /// Return an iterator the returns valid kmers
    pub fn kmers<'b, 'c>(&'b mut self, k: u8, canonical: bool) -> NuclKmer<'c>
    where
        'b: 'c,
    {
        if canonical {
            self.rev_seq = Some(self.seq.iter().rev().map(|n| complement(*n)).collect());
        }
        match self.rev_seq {
            Some(ref rev_seq) => NuclKmer::new(&self.seq, Some(&rev_seq), k),
            None => NuclKmer::new(&self.seq, None, k),
        }
    }

    /// Return an iterator the returns valid kmers in 4-bit form
    pub fn bit_kmers(&self, k: u8, canonical: bool) -> BitNuclKmer {
        BitNuclKmer::new(&self.seq, k, canonical)
    }

    /// Construct an owned version of `self` to, e.g. pass across threads
    /// (it's not clear why this can't be the `impl for Clone`, but the
    /// 'static lifetime doesn't work there for some reason)
    pub fn into_owned(self) -> SeqRecord<'static> {
        SeqRecord {
            id: Cow::Owned(self.id.clone().into_owned()),
            seq: Cow::Owned(self.seq.clone().into_owned()),
            qual: self.qual.clone().map(Cow::into_owned).map(Cow::Owned),
            rev_seq: self.rev_seq.clone(),
        }
    }
}

#[test]
fn test_quality_mask() {
    let seq_rec = SeqRecord {
        id: "".into(),
        // seq: Cow::Borrowed(&b"AGCT"[..]),
        seq: b"AGCT"[..].into(),
        qual: Some(b"AAA0"[..].into()),
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
            _ => unreachable!("Too many kmers"),
        }
        i += 1;
    }

    // test that we skip over N's
    i = 0;
    for (_, k, _) in SeqRecord::from_bytes(b"ACNGT").kmers(2, false) {
        match i {
            0 => assert_eq!(k, &b"AC"[..]),
            1 => assert_eq!(k, &b"GT"[..]),
            _ => unreachable!("Too many kmers"),
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
            _ => unreachable!("Too many kmers"),
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
            _ => unreachable!("Too many kmers"),
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
            _ => unreachable!("Too many kmers"),
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
            _ => unreachable!("Too many kmers"),
        }
        i += 1;
    }
}
