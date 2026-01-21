//! FAI (FASTA index) file format support for random access to FASTA files.
//!
//! The FAI format is a simple tab-delimited format that stores metadata about
//! sequences in a FASTA file, allowing random access to specific regions.
//!
//! Each line contains:
//! - NAME: sequence name
//! - LENGTH: length of the sequence in bases
//! - OFFSET: byte offset of the first base in the sequence
//! - LINEBASES: number of bases per line
//! - LINEWIDTH: number of bytes per line (including newline)
//!
//! # Example
//!
//! ```no_run
//! use needletail::parser::fai::FaiIndex;
//! use std::fs::File;
//! use std::io::{BufReader, Seek, SeekFrom, Read};
//!
//! let index = FaiIndex::from_file("sequence.fa.fai").expect("Failed to read index");
//!
//! // Get the region for chr1:100-200
//! if let Some(region) = index.fetch_region("chr1", 100, 200) {
//!     let mut file = File::open("sequence.fa").expect("Failed to open FASTA");
//!     file.seek(SeekFrom::Start(region.start_offset)).unwrap();
//!
//!     let mut buffer = vec![0u8; region.bytes_to_read];
//!     file.read_exact(&mut buffer).unwrap();
//!     // Process buffer, removing newlines
//! }
//! ```

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::errors::ParseError;

/// A single entry in the FAI index, representing one sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct FaiEntry {
    /// Name of the sequence
    pub name: String,
    /// Length of the sequence in bases
    pub length: u64,
    /// Byte offset of the first base in the sequence
    pub offset: u64,
    /// Number of bases per line (excluding newline)
    pub linebases: u64,
    /// Number of bytes per line (including newline)
    pub linewidth: u64,
}

impl FaiEntry {
    /// Parse a FAI entry from a tab-delimited line.
    ///
    /// # Errors
    ///
    /// Returns an error if the line doesn't have exactly 5 fields or if
    /// any of the numeric fields fail to parse.
    pub fn from_line(line: &str) -> Result<Self, ParseError> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            return Err(ParseError::new_invalid_fai_line(line.to_string()));
        }

        Ok(FaiEntry {
            name: fields[0].to_string(),
            length: fields[1]
                .parse()
                .map_err(|_| ParseError::new_invalid_fai_line(line.to_string()))?,
            offset: fields[2]
                .parse()
                .map_err(|_| ParseError::new_invalid_fai_line(line.to_string()))?,
            linebases: fields[3]
                .parse()
                .map_err(|_| ParseError::new_invalid_fai_line(line.to_string()))?,
            linewidth: fields[4]
                .parse()
                .map_err(|_| ParseError::new_invalid_fai_line(line.to_string()))?,
        })
    }

    /// Calculate the byte offset and number of bytes to read for a given region.
    ///
    /// # Arguments
    ///
    /// * `start` - 0-based start position in the sequence
    /// * `end` - 0-based end position (exclusive) in the sequence
    ///
    /// # Returns
    ///
    /// Returns `Some(FetchRegion)` if the region is valid, `None` if out of bounds.
    pub fn region(&self, start: u64, end: u64) -> Option<FetchRegion> {
        if start >= self.length || end > self.length || start >= end {
            return None;
        }

        // Calculate line and column for start position
        let start_line = start / self.linebases;
        let start_col = start % self.linebases;

        // Calculate line and column for end position
        let end_line = (end - 1) / self.linebases;
        let end_col = (end - 1) % self.linebases;

        // Calculate byte offsets
        let start_offset = self.offset + start_line * self.linewidth + start_col;
        let end_offset = self.offset + end_line * self.linewidth + end_col + 1;

        Some(FetchRegion {
            start_offset,
            bytes_to_read: (end_offset - start_offset) as usize,
            bases_to_fetch: (end - start) as usize,
            linebases: self.linebases,
            linewidth: self.linewidth,
        })
    }
}

/// Information needed to fetch a region from a FASTA file.
#[derive(Debug, Clone, PartialEq)]
pub struct FetchRegion {
    /// Byte offset in the FASTA file where reading should start
    pub start_offset: u64,
    /// Total number of bytes to read (may include newlines)
    pub bytes_to_read: usize,
    /// Number of bases that will be extracted (after removing newlines)
    pub bases_to_fetch: usize,
    /// Number of bases per line (for stripping newlines)
    pub linebases: u64,
    /// Number of bytes per line including newline (for stripping newlines)
    pub linewidth: u64,
}

impl FetchRegion {
    /// Remove newline characters from a buffer according to the FAI line format.
    ///
    /// This method filters out newline bytes and returns only the sequence bases.
    pub fn strip_newlines(&self, buffer: &[u8]) -> Vec<u8> {
        buffer
            .iter()
            .filter(|&&b| b != b'\n' && b != b'\r')
            .copied()
            .collect()
    }
}

/// A FAI index containing entries for all sequences in a FASTA file.
#[derive(Debug, Clone, Default)]
pub struct FaiIndex {
    /// Map from sequence name to index entry
    entries: HashMap<String, FaiEntry>,
    /// Sequence names in order
    names: Vec<String>,
}

impl FaiIndex {
    /// Create a new empty FAI index.
    pub fn new() -> Self {
        FaiIndex {
            entries: HashMap::new(),
            names: Vec::new(),
        }
    }

    /// Parse a FAI index from a reader.
    ///
    /// # Errors
    ///
    /// Returns an error if any line fails to parse.
    pub fn from_reader<R: Read>(reader: R) -> Result<Self, ParseError> {
        let buf_reader = BufReader::new(reader);
        let mut index = FaiIndex::new();

        for line_result in buf_reader.lines() {
            let line = line_result?;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }

            let entry = FaiEntry::from_line(trimmed)?;
            index.names.push(entry.name.clone());
            index.entries.insert(entry.name.clone(), entry);
        }

        Ok(index)
    }

    /// Parse a FAI index from a file path.
    ///
    /// # Errors
    ///
    /// Returns an error if the file cannot be opened or if parsing fails.
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, ParseError> {
        let file = File::open(path)?;
        Self::from_reader(file)
    }

    /// Get an entry by sequence name.
    pub fn get(&self, name: &str) -> Option<&FaiEntry> {
        self.entries.get(name)
    }

    /// Get the number of sequences in the index.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if the index is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Get an iterator over sequence names in order.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.names.iter().map(|s| s.as_str())
    }

    /// Get an iterator over all entries.
    pub fn entries(&self) -> impl Iterator<Item = &FaiEntry> {
        self.names.iter().filter_map(|name| self.entries.get(name))
    }

    /// Calculate the byte offset and number of bytes to read for a given region.
    ///
    /// # Arguments
    ///
    /// * `name` - sequence name
    /// * `start` - 0-based start position in the sequence
    /// * `end` - 0-based end position (exclusive) in the sequence
    ///
    /// # Returns
    ///
    /// Returns `Some(FetchRegion)` if the sequence exists and the region is valid,
    /// `None` otherwise.
    pub fn fetch_region(&self, name: &str, start: u64, end: u64) -> Option<FetchRegion> {
        self.entries.get(name).and_then(|e| e.region(start, end))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_fai_entry() {
        let line = "chr1\t248956422\t6\t70\t71";
        let entry = FaiEntry::from_line(line).unwrap();

        assert_eq!(entry.name, "chr1");
        assert_eq!(entry.length, 248956422);
        assert_eq!(entry.offset, 6);
        assert_eq!(entry.linebases, 70);
        assert_eq!(entry.linewidth, 71);
    }

    #[test]
    fn test_parse_fai_entry_invalid() {
        let line = "chr1\t248956422";
        let result = FaiEntry::from_line(line);
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_fai_index() {
        let fai_content = "chr1\t1000\t6\t70\t71\nchr2\t2000\t1085\t70\t71\n";
        let index = FaiIndex::from_reader(fai_content.as_bytes()).unwrap();

        assert_eq!(index.len(), 2);

        let chr1 = index.get("chr1").unwrap();
        assert_eq!(chr1.length, 1000);

        let chr2 = index.get("chr2").unwrap();
        assert_eq!(chr2.offset, 1085);
    }

    #[test]
    fn test_region_calculation() {
        // Sequence with 70 bases per line, 71 bytes per line (including \n)
        let entry = FaiEntry {
            name: "chr1".to_string(),
            length: 1000,
            offset: 6,
            linebases: 70,
            linewidth: 71,
        };

        // First 10 bases (within first line)
        let region = entry.region(0, 10).unwrap();
        assert_eq!(region.start_offset, 6);
        assert_eq!(region.bases_to_fetch, 10);

        // Bases 65-75 (spanning two lines)
        let region = entry.region(65, 75).unwrap();
        assert_eq!(region.bases_to_fetch, 10);
        // Should include the newline character
        assert!(region.bytes_to_read > 10);
    }

    #[test]
    fn test_region_out_of_bounds() {
        let entry = FaiEntry {
            name: "chr1".to_string(),
            length: 100,
            offset: 6,
            linebases: 70,
            linewidth: 71,
        };

        assert!(entry.region(0, 101).is_none()); // end beyond length
        assert!(entry.region(100, 101).is_none()); // start at length
        assert!(entry.region(50, 50).is_none()); // empty region
        assert!(entry.region(50, 40).is_none()); // start > end
    }

    #[test]
    fn test_strip_newlines() {
        let region = FetchRegion {
            start_offset: 0,
            bytes_to_read: 10,
            bases_to_fetch: 8,
            linebases: 4,
            linewidth: 5,
        };

        let buffer = b"ACGT\nTGCA";
        let stripped = region.strip_newlines(buffer);
        assert_eq!(stripped, b"ACGTTGCA");
    }

    #[test]
    fn test_index_names_order() {
        let fai_content = "seq3\t100\t0\t70\t71\nseq1\t200\t100\t70\t71\nseq2\t300\t300\t70\t71\n";
        let index = FaiIndex::from_reader(fai_content.as_bytes()).unwrap();

        let names: Vec<&str> = index.names().collect();
        assert_eq!(names, vec!["seq3", "seq1", "seq2"]);
    }
}
