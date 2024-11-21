use needletail::parse_fastx_file;

const TEST_FILES: [&str; 4] = [
    "./tests/data/test.fa.gz",
    "./tests/data/test.fa.bz2",
    "./tests/data/test.fa.xz",
    "./tests/data/test.fa.zst",
];

#[cfg(feature = "compression")]
#[test]
fn can_read_compressed_files_automatically() {
    use needletail::parser::Format;
    for p in &TEST_FILES {
        let mut reader = parse_fastx_file(p).unwrap();
        let mut i = 0;
        while let Some(record) = reader.next() {
            let seq = record.unwrap();
            assert_eq!(seq.format(), Format::Fasta);

            match i {
                0 => {
                    assert_eq!(seq.id(), b"test");
                    assert_eq!(seq.raw_seq(), b"AGCTGATCGA");
                    assert_eq!(seq.qual(), None);
                }
                1 => {
                    assert_eq!(seq.id(), b"test2");
                    assert_eq!(seq.raw_seq(), b"TAGC");
                    assert_eq!(seq.qual(), None);
                }
                _ => unreachable!("Too many records"),
            }
            i += 1;
        }
        assert_eq!(i, 2);
    }
}

#[cfg(not(feature = "compression"))]
#[test]
fn errors_on_compressed_files() {
    for p in &TEST_FILES {
        assert!(parse_fastx_file(p).is_err());
    }
}
