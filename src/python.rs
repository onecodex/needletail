//! Python bindings for needletail

// TODO:
// - Make the return values of `__repr__` and `__str__` show up as raw strings.
// - Make `normalize_seq`, `reverse_complement`, and `decode_phred` functions
//   able to handle `Record` objects as input.

use crate::parser::SequenceRecord;
use crate::quality::{decode_phred, PhredEncoding};
use crate::sequence::{complement, normalize};
use crate::{parse_fastx_file, parse_fastx_reader, FastxReader};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use pyo3::{create_exception, wrap_pyfunction};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::Cursor;
use std::path::PathBuf;

create_exception!(needletail, NeedletailError, pyo3::exceptions::PyException);

// Avoid some boilerplate with the error handling
macro_rules! py_try {
    ($call:expr) => {
        $call.map_err(|e| PyErr::new::<NeedletailError, _>(format!("{}", e)))?
    };
}

fn get_seq_snippet(seq: &str, max_len: usize) -> String {
    if seq.len() > max_len {
        let start = &seq[..max_len - 4];
        let end = &seq[seq.len() - 3..];
        format!("{}…{}", start, end)
    } else {
        seq.to_string()
    }
}

/// An iterator that yields sequence records.
///
/// Yields
/// ------
/// Record
///     A `Record` object representing a sequence record.
///
/// See also
/// --------
/// parse_fastx_file:
///     A function to parse sequence records from a FASTA/FASTQ file.
/// parse_fastx_string:
///     A function to parse sequence records from a FASTA/FASTQ string.
/// Record:
///     A class representing a FASTA/FASTQ sequence record.
#[pyclass]
pub struct PyFastxReader {
    reader: Box<dyn FastxReader>,
}

#[pymethods]
impl PyFastxReader {
    fn __repr__(&self) -> PyResult<String> {
        Ok("<FastxReader>".to_string())
    }

    fn __iter__(slf: PyRefMut<Self>) -> PyRefMut<Self> {
        slf
    }

    fn __next__(mut slf: PyRefMut<Self>) -> PyResult<Option<Record>> {
        if let Some(rec) = slf.reader.next() {
            let record = py_try!(rec);
            Ok(Some(Record::from_sequence_record(&record)))
        } else {
            Ok(None)
        }
    }
}

/// A record representing a biological sequence.
///
/// Parameters
/// ----------
/// id : str
///     The identifier of the sequence record.
/// seq : str
///     A string representing the sequence.
///
/// Attributes
/// ----------
/// id : str
///     The identifier of the sequence record. In a FASTA file, this is the
///     string containing all characters (including whitespaces) after the
///     leading '>' character. In a FASTQ file, this is the string containing
///     all characters (including whitespaces) after the leading '@' character.
/// seq : str
///     A string representing the sequence.
/// qual : str, optional
///     A string representing the quality scores of the sequence. If the object
///     represents a FASTA record, this attribute will be `None`.
/// name : str
///     The name of the sequence record. This is the string before the first
///     whitespace character in the `id` attribute.
/// description : str, optional
///     The description of the sequence record. This is the string after the
///     first whitespace character in the `id` attribute. If the `id` attribute
///     contains no whitespace characters, this attribute will be `None`.
///
/// Methods
/// -------
/// is_fasta
///     Check if the object represents a FASTA record.
/// is_fastq
///     Check if the object represents a FASTQ record.
/// normalize(iupac)
///     Normalize the sequence stored in the `seq` attribute of the object.
#[pyclass]
pub struct Record {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    seq: String,
    #[pyo3(get)]
    qual: Option<String>,
}

impl Record {
    fn from_sequence_record(rec: &SequenceRecord) -> Self {
        Self {
            id: String::from_utf8_lossy(rec.id()).to_string(),
            seq: String::from_utf8_lossy(&rec.seq()).to_string(),
            qual: rec.qual().map(|q| String::from_utf8_lossy(q).to_string()),
        }
    }
}

#[pymethods]
impl Record {
    #[getter]
    pub fn name(&self) -> PyResult<&str> {
        if let Some(pos) = self.id.find(char::is_whitespace) {
            Ok(&self.id[..pos])
        } else {
            Ok(&self.id)
        }
    }

    #[getter]
    pub fn description(&self) -> PyResult<Option<&str>> {
        if let Some(pos) = self.id.find(char::is_whitespace) {
            Ok(Some(&self.id[pos..].trim_start()))
        } else {
            Ok(None)
        }
    }

    /// Check if the object represents a FASTA record.
    ///
    /// Returns
    /// -------
    /// bool
    ///     `True` if the record lacks quality information, otherwise `False`.
    pub fn is_fasta(&self) -> PyResult<bool> {
        Ok(self.qual.is_none())
    }

    /// Check if the object represents a FASTQ record.
    ///
    /// Returns
    /// -------
    /// bool
    ///     `True` if the record has quality information, otherwise `False`.
    pub fn is_fastq(&self) -> PyResult<bool> {
        Ok(self.qual.is_some())
    }

    /// Normalize the sequence stored in the `seq` attribute of the object.
    ///
    /// See also
    /// --------
    /// normalize_seq: A function to normalize nucleotide sequence strings.
    #[pyo3(signature = (iupac=false))]
    pub fn normalize(&mut self, iupac: bool) -> PyResult<()> {
        if let Some(s) = normalize(self.seq.as_bytes(), iupac) {
            self.seq = String::from_utf8_lossy(&s).to_string();
        }
        Ok(())
    }

    #[new]
    #[pyo3(signature = (id, seq, qual=None))]
    fn new(id: String, seq: String, qual: Option<String>) -> PyResult<Record> {
        // If `qual` is not None, check if it has the same length as `seq`
        if let Some(qual) = &qual {
            if qual.len() != seq.len() {
                return Err(PyValueError::new_err(
                    "Sequence and quality strings must have the same length",
                ));
            }
        }
        Ok(Record { id, seq, qual })
    }

    pub fn __hash__(&self) -> PyResult<u64> {
        let mut hasher = DefaultHasher::new();
        self.id.hash(&mut hasher);
        self.seq.hash(&mut hasher);
        match &self.qual {
            Some(qual) => qual.hash(&mut hasher),
            None => {}
        }
        Ok(hasher.finish())
    }

    pub fn __eq__(&self, other: &Record) -> PyResult<bool> {
        Ok(self.id == other.id && self.seq == other.seq && self.qual == other.qual)
    }

    pub fn __len__(&self) -> PyResult<usize> {
        Ok(self.seq.len())
    }

    pub fn __str__(&self) -> PyResult<String> {
        if self.qual.is_none() {
            Ok(format!(">{}\n{}\n", self.id, self.seq))
        } else {
            Ok(format!(
                "@{}\n{}\n+\n{}\n",
                self.id,
                self.seq,
                self.qual.clone().unwrap()
            ))
        }
    }

    fn __repr__(&self) -> PyResult<String> {
        let id_snippet = match self.name() {
            Ok(name) if name != self.id => format!("{}…", name),
            Ok(name) => name.to_string(),
            Err(_) => self.id.clone(),
        };
        let seq_snippet = get_seq_snippet(&self.seq, 20);
        let quality_snippet = match &self.qual {
            Some(qual) => get_seq_snippet(qual, 20),
            None => "None".to_string(),
        };
        Ok(format!(
            "Record(id={}, seq={}, qual={})",
            id_snippet, seq_snippet, quality_snippet
        ))
    }
}

/// An iterator that reads sequence records from a FASTA/FASTQ file.
///
/// Parameters
/// ----------
/// path : str or pathlib.Path
///     The path to a FASTA/FASTQ file.
///
/// Returns
/// -------
/// PyFastxReader
///     A `PyFastxReader` iterator that yields `Record` objects representing
///     sequences from the input file.
///
/// Raises
/// ------
/// NeedletailError
///     If an error occurs while reading and parsing the input file.
///
/// See also
/// --------
/// parse_fastx_string:
///     A function to parse sequence records from a FASTA/FASTQ string.
/// PyFastxReader:
///     A class with instances that are iterators that yield `Record` objects.
#[pyfunction]
#[pyo3(name = "parse_fastx_file")]
fn py_parse_fastx_file(path: PathBuf) -> PyResult<PyFastxReader> {
    let reader = py_try!(parse_fastx_file(path));
    Ok(PyFastxReader { reader })
}

/// Parse sequence records from a FASTA/FASTQ string.
///
/// Parameters
/// ----------
/// content : str
///     A string containing FASTA/FASTQ-formatted sequence records.
///
/// Returns
/// -------
/// PyFastxReader
///     A `PyFastxReader` iterator that yields `Record` objects representing
///     sequences from the input string.
///
/// Raises
/// ------
/// NeedletailError
///     If an error occurs while parsing the input string.
///
/// See also
/// --------
/// parse_fastx_file:
///     A function to parse sequence records from a FASTA/FASTQ file.
/// PyFastxReader:
///     A class with instances that are iterators that yield `Record` objects.
#[pyfunction]
fn parse_fastx_string(content: &str) -> PyResult<PyFastxReader> {
    let reader = py_try!(parse_fastx_reader(Cursor::new(content.to_owned())));
    Ok(PyFastxReader { reader })
}

/// Normalize the sequence string of nucleotide records by:
///
/// - Converting lowercase characters to uppercase.
/// - Removing whitespace and newline characters.
/// - Replacing 'U' with 'T'.
/// - Replacing '.' and '~' with '-'.
/// - Replacing characters not in 'ACGTN-' with 'N', unless `iupac` is `True`,
///   in which case characters representing nucleotide ambiguity are not
///   replaced.
///
/// Parameters
/// ----------
/// seq : str
///     A string representing a nucleotide sequence.
/// iupac : bool, default: False
///     If `True`, characters representing nucleotide ambiguity ('B', 'D',
///     'H', 'V', 'R', 'Y', 'S', 'W', 'K', and 'M', and their lowercase
///     forms) will not be converted to 'N'. Lowercase characters will still
///     be converted to uppercase.
///
/// Returns
/// -------
/// str
///     The normalized sequence string.
///
/// Notes
/// -----
/// The `normalize` method is designed for nucleotide sequences only. If
/// used with protein sequences, it will incorrectly process amino acid
/// characters as if they were nucleotides.
#[pyfunction]
#[pyo3(signature = (seq, iupac=false))]
pub fn normalize_seq(seq: &str, iupac: bool) -> PyResult<String> {
    if let Some(s) = normalize(seq.as_bytes(), iupac) {
        Ok(String::from_utf8_lossy(&s).to_string())
    } else {
        Ok(seq.to_owned())
    }
}

/// Compute the reverse complement of a nucleotide sequence.
///
/// Parameters:
/// -----------
/// seq : str
///     A string representing a nucleotide sequence.
///
/// Returns:
/// --------
/// str
///     The reverse complement of the input nucleotide sequence.
#[pyfunction]
pub fn reverse_complement(seq: &str) -> PyResult<String> {
    let comp: Vec<u8> = seq
        .as_bytes()
        .iter()
        .rev()
        .map(|n| complement(*n))
        .collect();
    Ok(String::from_utf8_lossy(&comp).to_string())
}

/// Decode Phred quality data to quality scores.
///
/// Parameters:
/// -----------
/// phred : str
///     A string representing Phred-encoded quality data.
/// base_64 : bool, default=False
///     If `True`, return the quality using the Phred+64 encoding, otherwise
///     the Phred+33 encoding will be used.
///
/// Returns:
/// --------
/// tuple of int
///     A list of integers representing quality scores derived from the
///     probability of a base-calling error using a logarithmic transformation.
#[pyfunction]
#[pyo3(name = "decode_phred", signature = (qual, base_64=false))]
pub fn py_decode_phred(qual: &str, base_64: bool, py: Python<'_>) -> PyResult<Py<PyTuple>> {
    let encoding = if base_64 {
        PhredEncoding::Phred64
    } else {
        PhredEncoding::Phred33
    };
    decode_phred(qual.as_bytes(), encoding)
        .map(|vec| PyTuple::new_bound(py, &vec).into())
        .map_err(|e| PyValueError::new_err(format!("Invalid Phred quality: {}", e)))
}

#[pymodule]
fn needletail(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFastxReader>()?;
    m.add_class::<Record>()?;
    m.add_wrapped(wrap_pyfunction!(py_parse_fastx_file))?;
    m.add_wrapped(wrap_pyfunction!(parse_fastx_string))?;
    m.add_wrapped(wrap_pyfunction!(normalize_seq))?;
    m.add_wrapped(wrap_pyfunction!(reverse_complement))?;
    m.add_wrapped(wrap_pyfunction!(py_decode_phred))?;
    m.add("NeedletailError", py.get_type_bound::<NeedletailError>())?;
    Ok(())
}
