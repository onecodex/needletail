use std::fs;

use needletail::errors::ParseError;
use needletail::parser::parse_fastx_file;
use serde_derive::Deserialize;

#[derive(Debug, Deserialize)]
struct TestCase {
    filename: String,
    // origin: String,
    tags: Option<Vec<String>>,
    // comments: Option<Vec<String>>,
}

#[derive(Debug, Deserialize)]
struct TestIndex {
    valid: Vec<TestCase>,
    invalid: Option<Vec<TestCase>>,
}

fn test_fastx_file(path: &str) -> Result<(), ParseError> {
    let mut reader = parse_fastx_file(path)?;
    while let Some(rec) = reader.next() {
        let _ = rec?;
    }
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

        let path = format!("tests/specimen/FASTA/{}", test.filename);
        assert_eq!(test_fastx_file(&path), Ok(()));
    }
}

#[test]
fn test_specimen_fastq() {
    let raw_index = fs::read_to_string("tests/specimen/FASTQ/index.toml").unwrap();
    let index: TestIndex = toml::from_str(&raw_index).expect("Could not deserialize index");

    for test in index.valid {
        if test.filename == "wrapping_original_sanger.fastq"
            || test.filename == "longreads_original_sanger.fastq"
            || test.filename == "tricky.fastq"
        {
            // may god have mercy upon us if someone ever tries a file like this
            // (sequences are one-line, but quality scores are line-wrapped)
            continue;
        }

        let path = format!("tests/specimen/FASTQ/{}", test.filename);
        assert!(
            test_fastx_file(&path).is_ok(),
            "File {} is bad?",
            test.filename
        );
    }

    for test in index.invalid.unwrap_or_default() {
        if test.filename == "error_diff_ids.fastq" {
            // we don't care if the sequence ID doesn't match the quality id?
            continue;
        }

        // We don't check for ascii validity since it's a big hit perf wise
        // This means some invalid sequences are considered ok but it's not a big issue
        // in practice
        if test.filename.starts_with("error_qual_")
            || test.filename == "error_spaces.fastq"
            || test.filename == "error_tabs.fastq"
        {
            continue;
        }

        let path = format!("tests/specimen/FASTQ/{}", test.filename);
        assert!(
            test_fastx_file(&path).is_err(),
            "File {} is good?",
            test.filename
        );
    }
}
