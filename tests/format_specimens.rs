use std::fs;
use std::io::Read;

use needletail::formats::{FastaParser, FastqParser, RecParser};
use needletail::{ParseError, ParseErrorType};
use serde_derive::Deserialize;
use toml;

#[derive(Deserialize)]
struct TestCase {
    filename: String,
    // origin: String,
    tags: Option<Vec<String>>,
    // comments: Option<Vec<String>>,
}

#[derive(Deserialize)]
struct TestIndex {
    valid: Vec<TestCase>,
    invalid: Option<Vec<TestCase>>,
}

fn test_fasta_file(reader: &mut dyn Read, filename: &str) -> Result<(), ParseError> {
    let mut data: Vec<u8> = Vec::new();
    let _ = reader.read_to_end(&mut data)?;

    let parser = FastaParser::new(&data, true)
        .unwrap_or_else(|_| panic!("Can not open test data: {}", filename));
    let record_number = 0;
    for record in parser {
        let _ = record.map_err(|e| e.record(record_number))?;
    }
    Ok(())
}

fn test_fastq_file(reader: &mut dyn Read, filename: &str) -> Result<(), ParseError> {
    let mut data: Vec<u8> = Vec::new();
    let _ = reader.read_to_end(&mut data)?;

    let mut parser = FastqParser::new(&data, true)
        .unwrap_or_else(|_| panic!("Can not open test data: {}", filename));
    let record_number = 0;
    for record in parser.by_ref() {
        let rec = record.map_err(|e| e.record(record_number))?;
        if !rec
            .qual
            .iter()
            .all(|c| (*c >= b'!' && *c <= b'~') || *c == b'\n')
        {
            return Err(ParseError::new(
                "FASTQ has bad quality scores",
                ParseErrorType::Invalid,
            ));
        }
    }
    parser.eof()?;
    Ok(())
}

#[test]
fn test_specimen_fasta() {
    let raw_index = fs::read_to_string("tests/specimen/FASTA/index.toml").unwrap();
    let index: TestIndex = toml::from_str(&raw_index).expect("Could not deserialize index");
    for test in index.valid {
        // what kind of sicko puts comments in FASTAs?
        if test
            .tags
            .unwrap_or_else(Vec::new)
            .contains(&String::from("comments"))
        {
            continue;
        }

        let mut test_content =
            fs::File::open(&format!("tests/specimen/FASTA/{}", test.filename)).unwrap();
        assert_eq!(test_fasta_file(&mut test_content, &test.filename), Ok(()));
    }
}

#[test]
fn test_specimen_fastq() {
    let raw_index = fs::read_to_string("tests/specimen/FASTQ/index.toml").unwrap();
    let index: TestIndex = toml::from_str(&raw_index).expect("Could not deserialize index");

    for test in index.valid {
        if test.filename == "wrapping_original_sanger.fastq" {
            // may god have mercy upon us if someone ever tries a file like this
            // (sequences are one-line, but quality scores are line-wrapped)
            continue;
        }
        let mut test_content =
            fs::File::open(&format!("tests/specimen/FASTQ/{}", test.filename)).unwrap();
        assert_eq!(
            test_fastq_file(&mut test_content, &test.filename),
            Ok(()),
            "File {} is bad?",
            test.filename
        );
    }

    for test in index.invalid.unwrap_or_else(Vec::new) {
        if test.filename == "error_diff_ids.fastq" {
            // we don't care if the sequence ID doesn't match the quality id?
            continue;
        }
        let mut test_content =
            fs::File::open(&format!("tests/specimen/FASTQ/{}", test.filename)).unwrap();
        assert!(
            test_fastq_file(&mut test_content, &test.filename).is_err(),
            format!("File {} is good?", test.filename)
        );
    }
}
