//! For working with sequences that have identifiers and optionally quality
//! information.
//!
//! Primarily used as a common intermediate for processing both FASTA and
//! FASTQ data.
use std::borrow::Cow;
use std::io::Write;

use memchr::memchr;

use crate::sequence::{QualitySequence, Sequence};
use crate::util::ParseError;

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

/// An intermediate structure for handling sequence data and harmonizing both
/// FASTA and FASTQ records into a common format.
pub struct SequenceRecord<'a> {
    pub id: Cow<'a, [u8]>,
    pub seq: Cow<'a, [u8]>,
    pub qual: Option<Cow<'a, [u8]>>,
}

impl<'a> SequenceRecord<'a> {
    /// Creates a new SequenceRecord
    pub fn new(id: Cow<'a, [u8]>, seq: Cow<'a, [u8]>, qual: Option<Cow<'a, [u8]>>) -> Self {
        SequenceRecord { id, seq, qual }
    }

    /// Fixes up potential problems with sequence headers including tabs being
    /// present (may break downstream analyses with headers in TSVs) and with
    /// non-UTF8 characters being present, e.g. non-breaking spaces on Windows
    /// encodings (0x0A) breaks some tools.
    pub fn mask_header(mut self) -> Self {
        if let Some(id) = mask_header_tabs(&self.id) {
            self.id = id.into();
        }
        if let Some(id) = mask_header_utf8(&self.id) {
            self.id = id.into();
        }
        self
    }

    /// Write this SequenceRecord to writer as a FASTA with the provided line
    /// ending (ending should be either `\r\n` or preferably `\n`).
    pub fn write_fasta(&self, writer: &mut dyn Write, ending: &[u8]) -> Result<(), ParseError> {
        writer.write_all(b">")?;
        writer.write_all(&self.id)?;
        writer.write_all(ending)?;
        writer.write_all(&self.seq)?;
        writer.write_all(ending)?;
        Ok(())
    }

    /// Write this SequenceRecord to writer as a FASTQ with the provided line
    /// ending (ending should be either `\r\n` or preferably `\n`).
    pub fn write_fastq(&self, writer: &mut dyn Write, ending: &[u8]) -> Result<(), ParseError> {
        writer.write_all(b"@")?;
        writer.write_all(&self.id)?;
        writer.write_all(ending)?;
        writer.write_all(&self.seq)?;
        writer.write_all(ending)?;
        writer.write_all(b"+")?;
        writer.write_all(ending)?;
        // this is kind of a hack, but we want to allow writing out sequences
        // that don't have qualitys so this will mask to "good" if the quality
        // slice is empty
        if let Some(qual) = &self.qual {
            writer.write_all(&qual)?;
        } else {
            writer.write_all(&vec![b'I'; self.seq.len()])?;
        }
        writer.write_all(ending)?;
        Ok(())
    }
}

impl<'a> From<&'a [u8]> for SequenceRecord<'a> {
    fn from(slice: &'a [u8]) -> Self {
        SequenceRecord::new(Cow::from(&b""[..]), slice.into(), None)
    }
}

impl<'a> Sequence<'a> for SequenceRecord<'a> {
    fn sequence(&'a self) -> &'a [u8] {
        self.seq.as_ref()
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
