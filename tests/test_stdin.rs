use std::io::{Seek, SeekFrom, Write};

use assert_cmd::prelude::*;
use predicates::str::contains;

#[cfg(feature = "compression")]
#[test]
fn test_stdin_gz() {
    // Generated with `echo ">id1\nAGTCGTCA" | gzip -c | xxd  -i`
    let input: &[u8] = &[
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xb3, 0xcb, 0x4c, 0x31, 0xe4,
        0x72, 0x74, 0x0f, 0x71, 0x06, 0x22, 0x47, 0x2e, 0x00, 0x9e, 0x3a, 0x32, 0x8c, 0x0e, 0x00,
        0x00, 0x00,
    ];
    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(input).unwrap();
    file.flush().unwrap();
    file.seek(SeekFrom::Start(0)).unwrap();

    escargot::CargoBuild::new()
        .example("stdin_pipe")
        .current_release()
        .current_target()
        .run()
        .unwrap()
        .command()
        .stdin(file.into_file())
        .assert()
        .success()
        .stdout(contains("There are 8 bases in your file"))
        .stdout(contains("There are 0 AAAAs in your file"));
}

#[cfg(feature = "compression")]
#[test]
fn test_stdin_xz() {
    // Generated with `echo ">id1\nAGTCGTCA" | xz -c | xxd  -i`
    let input: &[u8] = &[
        0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00, 0x00, 0x04, 0xe6, 0xd6, 0xb4, 0x46, 0x02, 0x00, 0x21,
        0x01, 0x16, 0x00, 0x00, 0x00, 0x74, 0x2f, 0xe5, 0xa3, 0x01, 0x00, 0x0d, 0x3e, 0x69, 0x64,
        0x31, 0x0a, 0x41, 0x47, 0x54, 0x43, 0x47, 0x54, 0x43, 0x41, 0x0a, 0x00, 0x00, 0x00, 0x12,
        0x0f, 0x91, 0x75, 0xef, 0x7b, 0x63, 0x17, 0x00, 0x01, 0x26, 0x0e, 0x08, 0x1b, 0xe0, 0x04,
        0x1f, 0xb6, 0xf3, 0x7d, 0x01, 0x00, 0x00, 0x00, 0x00, 0x04, 0x59, 0x5a,
    ];
    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(input).unwrap();
    file.flush().unwrap();
    file.seek(SeekFrom::Start(0)).unwrap();

    escargot::CargoBuild::new()
        .example("stdin_pipe")
        .current_release()
        .current_target()
        .run()
        .unwrap()
        .command()
        .stdin(file.into_file())
        .assert()
        .success()
        .stdout(contains("There are 8 bases in your file"))
        .stdout(contains("There are 0 AAAAs in your file"));
}

#[cfg(feature = "compression")]
#[test]
fn test_stdin_bzip() {
    // Generated with `echo ">id1\nAGTCGTCA" | bzip2 -c | xxd  -i`
    let input: &[u8] = &[
        0x42, 0x5a, 0x68, 0x39, 0x31, 0x41, 0x59, 0x26, 0x53, 0x59, 0x9f, 0x9d, 0xf9, 0xa2, 0x00,
        0x00, 0x01, 0xcf, 0x00, 0x00, 0x10, 0x20, 0x01, 0x28, 0x80, 0x04, 0x00, 0x04, 0x20, 0x20,
        0x00, 0x22, 0x0c, 0x9a, 0x64, 0x20, 0xc9, 0x88, 0x21, 0x95, 0x8e, 0x82, 0x75, 0x27, 0x8b,
        0xb9, 0x22, 0x9c, 0x28, 0x48, 0x4f, 0xce, 0xfc, 0xd1, 0x00,
    ];
    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(input).unwrap();
    file.flush().unwrap();
    file.seek(SeekFrom::Start(0)).unwrap();

    escargot::CargoBuild::new()
        .example("stdin_pipe")
        .current_release()
        .current_target()
        .run()
        .unwrap()
        .command()
        .stdin(file.into_file())
        .assert()
        .success()
        .stdout(contains("There are 8 bases in your file"))
        .stdout(contains("There are 0 AAAAs in your file"));
}

#[cfg(feature = "compression")]
#[test]
fn test_stdin_zstd() {
    // Generated with `echo ">id1\nAGTCGTCA" | zstd -c | xxd  -i`
    let input: &[u8] = &[
        0x28, 0xb5, 0x2f, 0xfd, 0x04, 0x58, 0x71, 0x00, 0x00, 0x3e, 0x69, 0x64, 0x31, 0x0a, 0x41,
        0x47, 0x54, 0x43, 0x47, 0x54, 0x43, 0x41, 0x0a, 0x52, 0x9d, 0x37, 0x8d,
    ];
    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(input).unwrap();
    file.flush().unwrap();
    file.seek(SeekFrom::Start(0)).unwrap();

    escargot::CargoBuild::new()
        .example("stdin_pipe")
        .current_release()
        .current_target()
        .run()
        .unwrap()
        .command()
        .stdin(file.into_file())
        .assert()
        .success()
        .stdout(contains("There are 8 bases in your file"))
        .stdout(contains("There are 0 AAAAs in your file"));
}

#[test]
fn test_stdin_no_compression() {
    let input: &[u8] = b">id1\nAGTCGTCA";
    let mut file = tempfile::NamedTempFile::new().unwrap();
    file.write_all(input).unwrap();
    file.flush().unwrap();
    file.seek(SeekFrom::Start(0)).unwrap();

    escargot::CargoBuild::new()
        .example("stdin_pipe")
        .current_release()
        .current_target()
        .run()
        .unwrap()
        .command()
        .stdin(file.into_file())
        .assert()
        .success()
        .stdout(contains("There are 8 bases in your file"))
        .stdout(contains("There are 0 AAAAs in your file"));
}
