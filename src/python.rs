//! Python bindings for needletail

use std::io::Cursor;

use pyo3::prelude::*;
use pyo3::{create_exception, wrap_pyfunction};
use pyo3::{PyIterProtocol, PyObjectProtocol};

use crate::sequence::{complement, normalize};
use crate::{
    parse_fastx_file as rs_parse_fastx_file, parse_fastx_reader, parser::SequenceRecord,
    FastxReader,
};

create_exception!(needletail, NeedletailError, pyo3::exceptions::Exception);

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

#[pyproto]
impl PyObjectProtocol for PyFastxReader {
    fn __repr__(&self) -> PyResult<String> {
        Ok("<FastxParser>".to_string())
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

#[pyclass]
pub struct FastxReaderIterator {
    t: PyObject,
}

#[pyproto]
impl PyIterProtocol for PyFastxReader {
    fn __iter__(slf: PyRefMut<Self>) -> PyResult<FastxReaderIterator> {
        let gil_guard = Python::acquire_gil();
        let py = gil_guard.python();
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
    pub fn is_fasta(&self) -> PyResult<bool> {
        Ok(self.qual.is_none())
    }

    pub fn is_fastq(&self) -> PyResult<bool> {
        Ok(self.qual.is_some())
    }

    pub fn normalize(&mut self, iupac: bool) -> PyResult<()> {
        if let Some(s) = normalize(self.seq.as_bytes(), iupac) {
            self.seq = String::from_utf8_lossy(&s).to_string();
        }
        Ok(())
    }
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
pub fn reverse_complement(seq: &str) -> PyResult<String> {
    let comp: Vec<u8> = seq
        .as_bytes()
        .iter()
        .rev()
        .map(|n| complement(*n))
        .collect();
    Ok(String::from_utf8_lossy(&comp).to_string())
}

#[pyproto]
impl PyIterProtocol for FastxReaderIterator {
    fn __next__(slf: PyRef<Self>) -> PyResult<Option<Record>> {
        let gil_guard = Python::acquire_gil();
        let py = gil_guard.python();
        let mut parser: PyRefMut<PyFastxReader> = slf.t.extract(py)?;
        if let Some(rec) = parser.reader.next() {
            let record = py_try!(rec);
            Ok(Some(Record::from_sequence_record(&record)))
        } else {
            Ok(None)
        }
    }
}

#[pymodule]
fn needletail(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyFastxReader>()?;
    m.add_wrapped(wrap_pyfunction!(parse_fastx_file))?;
    m.add_wrapped(wrap_pyfunction!(parse_fastx_string))?;
    m.add_wrapped(wrap_pyfunction!(normalize_seq))?;
    m.add_wrapped(wrap_pyfunction!(reverse_complement))?;
    m.add("NeedletailError", py.get_type::<NeedletailError>())?;

    Ok(())
}
