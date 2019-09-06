use std::io::Read;

use needletail::formats::{FastaParser, FastqParser, RecParser};
use needletail::{ParseError, ParseErrorType};
use reqwest::get;
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

#[test]
#[ignore]
fn test_specimen_fasta() {
    let base_path = "https://raw.githubusercontent.com/BioJulia/FormatSpecimens.jl/master/FASTA";
    let idx_path = format!("{}/index.toml", base_path);
    let raw_index = get(&idx_path)
        .expect("Could not retrieve index")
        .text()
        .expect("Could not decode index");

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

        let test_path = format!("{}/{}", base_path, test.filename);
        let mut test_reader = get(&test_path).expect("Could not retrieve test data");
        assert_eq!(test_fasta_file(&mut test_reader, &test.filename), Ok(()));
    }
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
#[ignore]
fn test_specimen_fastq() {
    let base_path = "https://raw.githubusercontent.com/BioJulia/FormatSpecimens.jl/master/FASTQ/";
    let idx_path = format!("{}/index.toml", base_path);
    let raw_index = get(&idx_path)
        .expect("Could not retrieve index")
        .text()
        .expect("Could not decode index");

    let index: TestIndex = toml::from_str(&raw_index).expect("Could not deserialize index");

    for test in index.valid {
        if test.filename == "wrapping_original_sanger.fastq" {
            // may god have mercy upon us if someone ever tries a file like this
            // (sequences are one-line, but quality scores are line-wrapped)
            continue;
        }
        let test_path = format!("{}/{}", base_path, test.filename);
        let mut test_reader = get(&test_path).expect("Could not retrieve test data");
        assert_eq!(
            test_fastq_file(&mut test_reader, &test.filename),
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
        let test_path = format!("{}/{}", base_path, test.filename);
        let mut test_reader = get(&test_path).expect("Could not retrieve test data");
        assert!(
            test_fastq_file(&mut test_reader, &test.filename).is_err(),
            format!("File {} is good?", test.filename)
        );
    }
}
