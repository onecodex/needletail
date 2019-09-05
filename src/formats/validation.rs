use std::borrow::Cow;
use std::collections::HashSet;
use std::fmt;
use std::io::Write;
use std::iter::Iterator;

use crate::formats::{Fasta, Fastq, FastaParser, FastqParser, RecParser};
use crate::seq::{Sequence, mask_header_tabs, mask_header_utf8};
use crate::util::{ParseError, ParseErrorType};

#[derive(Debug, PartialEq, Eq, Hash, Copy, Clone)]
pub enum Coercion {
    /// Tab character in header
    TabInHeader,
    /// Bad character points in header
    NonUnicodeHeader,
    WhitespaceInSequence,
    GapInSequence,
    UridineInSequence,
    IupacInSequence,
    BadSequenceCharacter,
    TrimFastqSeparatorLine,
    /// An empty sequence has been removed
    EmptySequence,
    /// An empty newline was removed
    ExtraNewline,
}

impl fmt::Display for Coercion {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let msg = match &self {
            Coercion::TabInHeader => "Tab in header",
            Coercion::NonUnicodeHeader => "Header is not valid UTF8",
            Coercion::WhitespaceInSequence => "Whitespace in sequence",
            Coercion::GapInSequence => "Gaps in sequence",
            Coercion::UridineInSequence => "Uridine in sequence",
            Coercion::IupacInSequence => "IUPAC character in sequence",
            Coercion::BadSequenceCharacter => "Bad characters in sequence",
            Coercion::TrimFastqSeparatorLine => "FASTQ quality header was not empty",
            Coercion::EmptySequence => "Sequence was empty",
            Coercion::ExtraNewline => "Extra new lines removed",
        };
        write!(f, "{}", msg)
    }
}

impl From<Coercion> for ParseError {
    fn from(coercion: Coercion) -> ParseError {
        let err_type = match coercion {
            Coercion::TabInHeader | Coercion::NonUnicodeHeader => ParseErrorType::InvalidHeader,
            _ => ParseErrorType::InvalidRecord,
        };
        ParseError::new(format!("{}", coercion), err_type)
    }
}

fn validate_sequence(seq: &mut Sequence) -> (bool, Vec<Coercion>) {
    // TODO: we might want to add `context` to some of these errors?
    let mut errors = Vec::new();
    let mut skip_record = false;

    if let Some(s) = mask_header_utf8(&seq.id) {
        seq.id = Cow::Owned(s);
        errors.push(Coercion::NonUnicodeHeader);
    }

    if let Some(s) = mask_header_tabs(&seq.id) {
        seq.id = Cow::Owned(s);
        errors.push(Coercion::TabInHeader);
    }

    if seq.seq.len() == 0 {
        skip_record = true;
        errors.push(Coercion::EmptySequence);
    }

    let mut seq_changed = false;
    let new_seq = seq.seq.iter().filter_map(|n| {
        match *n {
            c @ b'A' | c @ b'a' | c @ b'C' | c @ b'c' | c @ b'G' | c @ b'g' | c @ b'T' | c @ b't' | c @ b'N' | c @ b'n' => Some(c),
            b' ' | b'\t' => {
                errors.push(Coercion::WhitespaceInSequence);
                seq_changed = true;
                None
            },
            b'-' | b'.' | b'~' => {
                errors.push(Coercion::GapInSequence);
                seq_changed = true;
                None
            },
            b'u' => {
                errors.push(Coercion::UridineInSequence);
                seq_changed = true;
                Some(b't')
            },
            b'U' => {
                errors.push(Coercion::UridineInSequence);
                seq_changed = true;
                Some(b'T')
            },

            b'B' | b'D' | b'H' | b'V' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' |
            b'b' | b'd' | b'h' | b'v' | b'r' | b'y' | b's' | b'w' | b'k' | b'm' => {
                errors.push(Coercion::IupacInSequence);
                seq_changed = true;
                Some(b'N')
            },
            _ => {
                errors.push(Coercion::BadSequenceCharacter);
                seq_changed = true;
                Some(b'N')
            }
        }
    }).collect();
    if seq_changed {
        seq.seq = Cow::Owned(new_seq);
    }

    (skip_record, errors)
}

pub struct FastqValidator {
    buffer: Vec<u8>,
    record: usize,
    strict: bool,
    coercions: HashSet<Coercion>,
}

impl FastqValidator {
    pub fn new(strict: bool) -> Result<Self, ParseError> {
        // we're not doing the full header handling here because we're not
        // dealing with formats that use headers (e.g. twobit) yet
        Ok(FastqValidator {
            buffer: Vec::new(),
            record: 0,
            strict,
            coercions: HashSet::new(),
        })
    }

    pub fn check_chunk<V>(
        &mut self,
        rec_writer: &mut dyn Write,
        data: &[u8],
        last: bool,
    ) -> Result<(), ParseError> {
        let mut all_coercions = HashSet::new();
        self.buffer.extend_from_slice(data);
        let used = {
            let mut rec_reader = FastqParser::from_buffer(&self.buffer, last);
            let mut start_pos = rec_reader.used();
            loop {
                let iter = rec_reader.by_ref();
                let write_raw = match iter.next() {
                    Some(Ok(s)) => {
                        let empty_quality_header = s.id2.is_empty();
                        let mut seq = Sequence::from(s);
                        // check sequence length
                        let (skip, mut coercions) = validate_sequence(&mut seq);
                        if !empty_quality_header {
                            // we automatically strip this out during 
                            // Fastq->Seq->Fastq conversion
                            coercions.push(Coercion::TrimFastqSeparatorLine);
                        }
                        if !coercions.is_empty() {
                            if self.strict {
                                // just return the first major issue?
                                return Err(coercions[0].into());
                            } else {
                                coercions.into_iter().for_each(|c| {
                                    all_coercions.insert(c);
                                });
                            }
                            if !skip {
                                Fastq::from(&seq).write(rec_writer)?;
                            } 
                            // TODO: add the coercions to the global coercions
                            false
                        } else {
                            true
                        }
                    },
                    Some(Err(e)) => Err(e.record(self.record))?,
                    None => break,
                };
                let end_pos = rec_reader.used();
                if write_raw {
                    rec_writer.write(&self.buffer[start_pos..end_pos])?;
                }
                start_pos = end_pos;
            }
            rec_reader.used()
        };
        self.buffer.drain(0..used);
        if last {
            let rec_reader = FastqParser::from_buffer(&self.buffer, last);
            rec_reader.eof().map_err(|e| e.record(self.record + 1))?;
        }
        all_coercions.into_iter().for_each(|c| {
            self.coercions.insert(c);
        });
        Ok(())
    }
}

pub struct FastaValidator {
    buffer: Vec<u8>,
    record: usize,
    strict: bool,
    coercions: HashSet<Coercion>,
}

impl FastaValidator {
    pub fn new(strict: bool) -> Result<Self, ParseError> {
        // we're not doing the full header handling here because we're not
        // dealing with formats that use headers (e.g. twobit) yet
        Ok(FastaValidator {
            buffer: Vec::new(),
            record: 0,
            strict,
            coercions: HashSet::new(),
        })
    }

    pub fn check_chunk<V>(
        &mut self,
        rec_writer: &mut dyn Write,
        data: &[u8],
        last: bool,
    ) -> Result<(), ParseError> {
        let mut all_coercions = HashSet::new();
        self.buffer.extend_from_slice(data);
        let used = {
            let mut rec_reader = FastaParser::from_buffer(&self.buffer, last);
            let mut start_pos = rec_reader.used();
            loop {
                let iter = rec_reader.by_ref();
                let write_raw = match iter.next() {
                    Some(Ok(s)) => {
                        let mut seq = Sequence::from(s);
                        // check sequence length
                        let (skip, coercions) = validate_sequence(&mut seq);
                        if !coercions.is_empty() {
                            if self.strict {
                                // just return the first major issue?
                                return Err(coercions[0].into());
                            } else {
                                coercions.into_iter().for_each(|c| {
                                    all_coercions.insert(c);
                                });
                            }
                            if !skip {
                                Fasta::from(&seq).write(rec_writer)?;
                            } 
                            false
                        } else {
                            true
                        }
                    },
                    Some(Err(e)) => Err(e.record(self.record))?,
                    None => break,
                };
                let end_pos = rec_reader.used();
                if write_raw {
                    rec_writer.write(&self.buffer[start_pos..end_pos])?;
                }
                start_pos = end_pos;
            }
            rec_reader.used()
        };
        self.buffer.drain(0..used);
        if last {
            let rec_reader = FastaParser::from_buffer(&self.buffer, last);
            rec_reader.eof().map_err(|e| e.record(self.record + 1))?;
        }
        all_coercions.into_iter().for_each(|c| {
            self.coercions.insert(c);
        });
        Ok(())
    }
}
