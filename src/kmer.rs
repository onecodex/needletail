//! Functions for splitting sequences into fixed-width moving windows (kmers)
//! and utilities for dealing with these kmers.

/// Returns true if the base is a unambiguous nucleic acid base (e.g. ACGT) and
/// false otherwise.
fn is_good_base(chr: u8) -> bool {
    match chr as char {
        'a' | 'c' | 'g' | 't' | 'A' | 'C' | 'G' | 'T' => true,
        _ => false,
    }
}

/// Generic moving window iterator over sequences to return k-mers
///
/// Iterator returns slices to the original data.
pub struct Kmers<'a> {
    k: u8,
    start_pos: usize,
    buffer: &'a [u8],
}

impl<'a> Kmers<'a> {
    /// Creates a new kmer-izer for a nucleotide/amino acid sequence.
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

/// A kmer-izer for a nucleotide acid sequences to return canonical kmers.
///
/// Iterator returns the position of the kmer, a slice to the original data,
/// and an boolean indicating if the kmer returned is the original or the
/// reverse complement.
pub struct CanonicalKmers<'a> {
    k: u8,
    start_pos: usize,
    buffer: &'a [u8],
    rc_buffer: &'a [u8],
}

impl<'a> CanonicalKmers<'a> {
    /// Creates a new iterator.
    ///
    /// It's generally more useful to use this directly from a sequences (e.g.
    /// `seq.canonical_kmers`. Requires a reference to the reverse complement
    /// of the sequence it's created on, e.g.
    /// ```
    /// use needletail::Sequence;
    /// use needletail::kmer::CanonicalKmers;
    ///
    /// let seq = b"ACGT";
    /// let rc = seq.reverse_complement();
    /// let c_iter = CanonicalKmers::new(seq, &rc, 3);
    /// for (pos, kmer, canonical) in c_iter {
    ///    // process data in here
    /// }
    ///
    /// ```
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

#[cfg(tests)]
mod tests {
    use super::*;

    #[test]
    fn can_kmerize() {
        let k_iter = Kmers::new(b"AGCT", 1);
        // test general function
        for (i, k) in k_iter.enumerate() {
            match i {
                0 => assert_eq!(k, &b"A"[..]),
                1 => assert_eq!(k, &b"G"[..]),
                2 => assert_eq!(k, &b"C"[..]),
                3 => assert_eq!(k, &b"T"[..]),
                _ => unreachable!("Too many kmers"),
            }
        }

        // test that we handle length 2 (and don't drop Ns)
        let k_iter = Kmers::new(b"AGNCT", 2);
        for (i, k) in k_iter.enumerate() {
            match i {
                0 => assert_eq!(k, &b"AC"[..]),
                1 => assert_eq!(k, &b"CN"[..]),
                2 => assert_eq!(k, &b"NG"[..]),
                3 => assert_eq!(k, &b"GT"[..]),
                _ => unreachable!("Too many kmers"),
            }
        }

        // test that the minimum length works
        let k_iter = Kmers::new(b"AC", 2);
        for k in k_iter {
            assert_eq!(k, &b"AC"[..]);
        }
    }

    #[test]
    fn can_canonicalize() {
        // test general function
        let seq = b"AGCT";
        let c_iter = CanonicalKmers(seq, &seq.reverse_complement(), 1);
        for (i, (_, k, is_c)) in c_iter.enumerate() {
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
        let c_iter = CanonicalKmers(seq, &seq.reverse_complement(), 2);
        for (i, (_, k, _)) in c_iter.enumerate() {
            match i {
                0 => assert_eq!(k, &b"AG"[..]),
                1 => assert_eq!(k, &b"GC"[..]),
                2 => assert_eq!(k, &b"AG"[..]),
                3 => assert_eq!(k, &b"TA"[..]),
                _ => unreachable!("Too many kmers"),
            }
        }

        let seq = b"AGNTA";
        let c_iter = CanonicalKmers(seq, &seq.reverse_complement(), 2);
        for (i, (ix, k, _)) in c_iter.enumerate() {
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
}
