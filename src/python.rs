//! Python bindings for needletail

use crate::sequence::{complement, normalize};
use crate::{
    parse_fastx_file as rs_parse_fastx_file, parse_fastx_reader, parser::SequenceRecord,
    FastxReader,
};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::{create_exception, wrap_pyfunction};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::io::Cursor;

create_exception!(needletail, NeedletailError, pyo3::exceptions::PyException);

// Avoid some boilerplate with the error handling
macro_rules! py_try {
    ($call:expr) => {
        $call.map_err(|e| PyErr::new::<NeedletailError, _>(format!("{}", e)))?
    };
}

#[pyclass]
pub struct PyFastxReader {
    reader: Box<dyn FastxReader>,
}

fn get_string_snippet(seq: &str, max_len: usize) -> String {
    if seq.len() > max_len {
        let start = &seq[..max_len - 4];
        let end = &seq[seq.len() - 3..];
        format!("{}…{}", start, end)
    } else {
        seq.to_string()
    }
}

#[pymethods]
impl PyFastxReader {
    fn __repr__(&self) -> PyResult<String> {
        Ok("<FastxParser>".to_string())
    }

    fn __iter__(slf: PyRefMut<Self>, py: Python<'_>) -> PyResult<FastxReaderIterator> {
        Ok(FastxReaderIterator { t: slf.into_py(py) })
    }
}

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
    pub fn is_fasta(&self) -> PyResult<bool> {
        Ok(self.qual.is_none())
    }

    #[getter]
    pub fn is_fastq(&self) -> PyResult<bool> {
        Ok(self.qual.is_some())
    }

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
            let wrapped_seq = self
                .seq
                .as_bytes()
                .chunks(60)
                .map(|chunk| String::from_utf8_lossy(chunk).to_string())
                .collect::<Vec<String>>()
                .join("\n");
            Ok(format!(">{}\n{}", self.id, wrapped_seq))
        } else {
            Ok(format!(
                "@{}\n{}\n+\n{}",
                self.id,
                self.seq,
                self.qual.clone().unwrap()
            ))
        }
    }

    fn __repr__(&self) -> PyResult<String> {
        let seq_snippet = get_string_snippet(&self.seq, 30);
        let quality_snippet = match &self.qual {
            Some(qual) => get_string_snippet(qual, 30),
            None => "None".to_string(),
        };
        Ok(format!(
            "Record(id={}, sequence={}, quality={})",
            self.id, seq_snippet, quality_snippet
        ))
    }
}

#[pyclass]
pub struct FastxReaderIterator {
    t: PyObject,
}

#[pymethods]
impl FastxReaderIterator {
    fn __next__(slf: PyRef<Self>, py: Python<'_>) -> PyResult<Option<Record>> {
        let mut parser: PyRefMut<PyFastxReader> = slf.t.extract(py)?;
        if let Some(rec) = parser.reader.next() {
            let record = py_try!(rec);
            Ok(Some(Record::from_sequence_record(&record)))
        } else {
            Ok(None)
        }
    }
}

// TODO: what would be really nice is to detect the type of pyobject so it would on file object etc
// not for initial release though

#[pyfunction]
fn parse_fastx_file(path: &str) -> PyResult<PyFastxReader> {
    let reader = py_try!(rs_parse_fastx_file(path));
    Ok(PyFastxReader { reader })
}

#[pyfunction]
fn parse_fastx_string(content: &str) -> PyResult<PyFastxReader> {
    let reader = py_try!(parse_fastx_reader(Cursor::new(content.to_owned())));
    Ok(PyFastxReader { reader })
}

#[pyfunction]
pub fn normalize_seq(seq: &str, iupac: bool) -> PyResult<String> {
    if let Some(s) = normalize(seq.as_bytes(), iupac) {
        Ok(String::from_utf8_lossy(&s).to_string())
    } else {
        Ok(seq.to_owned())
    }
}

#[pyfunction]
pub fn reverse_complement(seq: &str) -> String {
    let comp: Vec<u8> = seq
        .as_bytes()
        .iter()
        .rev()
        .map(|n| complement(*n))
        .collect();
    String::from_utf8_lossy(&comp).to_string()
}

#[pymodule]
fn needletail(py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyFastxReader>()?;
    m.add_class::<Record>()?;
    m.add_wrapped(wrap_pyfunction!(parse_fastx_file))?;
    m.add_wrapped(wrap_pyfunction!(parse_fastx_string))?;
    m.add_wrapped(wrap_pyfunction!(normalize_seq))?;
    m.add_wrapped(wrap_pyfunction!(reverse_complement))?;
    m.add("NeedletailError", py.get_type_bound::<NeedletailError>())?;
    Ok(())
}
