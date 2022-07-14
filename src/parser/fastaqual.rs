use std::borrow::Cow;

use crate::errors::ParseError;
use crate::parser::record::SequenceRecord;
use crate::parser::utils::{FastxReader, Format, LineEnding};

pub struct Reader<'a> {
    fasta_reader: Box<dyn FastxReader + 'a>,
    qual_reader: Box<dyn FastxReader + 'a>,
}

impl<'a> Reader<'a> {
    pub(crate) fn new(
        fasta_reader: Box<dyn FastxReader + 'a>,
        qual_reader: Box<dyn FastxReader + 'a>,
    ) -> Self {
        Self {
            fasta_reader,
            qual_reader,
        }
    }

    pub fn next(&mut self) -> Option<Result<FastaqualSequenceRecord, ParseError>> {
        match (self.fasta_reader.next(), self.qual_reader.next()) {
            (None, None) => None,
            (None, _) | (_, None) => Some(Err(ParseError::new_record_mismatch())),
            (Some(Err(e)), _) | (_, Some(Err(e))) => Some(Err(e)),
            (Some(Ok(fasta_record)), Some(Ok(qual_record))) => {
                match (fasta_record.format(), qual_record.format()) {
                    (Format::Fasta, Format::Fasta) => {
                        if fasta_record.id() == qual_record.id() {
                            Some(Ok(FastaqualSequenceRecord {
                                fasta_record,
                                qual_record,
                            }))
                        } else {
                            Some(Err(ParseError::new_record_mismatch()))
                        }
                    }
                    _ => Some(Err(ParseError::new_format_mismatch())),
                }
            }
        }
    }
}

#[derive(Debug)]
pub struct FastaqualSequenceRecord<'a> {
    fasta_record: SequenceRecord<'a>,
    qual_record: SequenceRecord<'a>,
}

impl<'a> FastaqualSequenceRecord<'a> {
    #[inline]
    pub fn format(&self) -> Format {
        Format::Fastaqual
    }

    #[inline]
    pub fn id(&self) -> &[u8] {
        self.fasta_record.id()
    }

    #[inline]
    pub fn raw_seq(&self) -> &[u8] {
        self.fasta_record.raw_seq()
    }

    pub fn seq(&self) -> Cow<[u8]> {
        self.fasta_record.seq()
    }

    #[inline]
    pub fn qual(&self) -> Result<Cow<[u8]>, ParseError> {
        let raw_qual = self.qual_record.raw_seq();
        let num_bases = self.num_bases();
        let mut new_buf = Vec::with_capacity(num_bases);
        let pieces = raw_qual.split(|e| e.is_ascii_whitespace());

        for piece in pieces {
            if piece.len() == 0 {
                continue;
            }

            let s = match std::str::from_utf8(piece) {
                Ok(s) => s,
                Err(_) => return Err(ParseError::new_invalid_qual_score()),
            };
            let score: u8 = match s.parse() {
                Ok(score) => score,
                Err(_) => return Err(ParseError::new_invalid_qual_score()),
            };
            new_buf.push(score);
        }

        if num_bases != new_buf.len() {
            return Err(ParseError::new_record_mismatch());
        }
        Ok(new_buf.into())
    }

    #[inline]
    pub fn all(&self) -> &[u8] {
        self.fasta_record.all()
    }

    #[inline]
    pub fn num_bases(&self) -> usize {
        self.fasta_record.num_bases()
    }

    pub fn start_line_number(&self) -> u64 {
        self.fasta_record.start_line_number()
    }

    pub fn line_ending(&self) -> LineEnding {
        self.fasta_record.line_ending()
    }
}
