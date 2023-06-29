use needletail::errors::IndexErrorKind;
use needletail::parser::IndexedReader;

#[test]
fn test_faidx() {
    let mut fai_reader = IndexedReader::from_path("tests/data/test_index.fa").unwrap();

    let subseq1 = fai_reader.subseq("xxx", Some(2), Some(8)).unwrap();
    assert_eq!(subseq1.seq, b"AAAGGG");
    assert_eq!(subseq1.name, "xxx");
    assert_eq!(subseq1.start, 2);
    assert_eq!(subseq1.end, 8);

    let subseq2 = fai_reader.subseq("xxx", Some(2), None).unwrap();
    assert_eq!(subseq2.seq, b"AAAGGGGG");
    assert_eq!(subseq2.name, "xxx");
    assert_eq!(subseq2.start, 2);
    assert_eq!(subseq2.end, 10);

    let subseq3 = fai_reader.subseq("xxx", None, Some(8)).unwrap();
    assert_eq!(subseq3.seq, b"AAAAAGGG");
    assert_eq!(subseq3.name, "xxx");
    assert_eq!(subseq3.start, 0);
    assert_eq!(subseq3.end, 8);

    let subseq4 = fai_reader.subseq("xxx", None, None).unwrap();
    assert_eq!(subseq4.seq, b"AAAAAGGGGG");
    assert_eq!(subseq4.name, "xxx");
    assert_eq!(subseq4.start, 0);
    assert_eq!(subseq4.end, 10);
    println!("{}", subseq4);

    let unknown_seq = fai_reader.subseq("unknown", None, None);
    assert_eq!(
        unknown_seq.err().unwrap().kind,
        IndexErrorKind::UnknownSeqName
    );

    let invalid_region1 = fai_reader.subseq("xxx", Some(8), Some(2));
    assert_eq!(
        invalid_region1.err().unwrap().kind,
        IndexErrorKind::InvalidRegion
    );

    let invalid_region2 = fai_reader.subseq("xxx", Some(2), Some(6666));
    assert_eq!(
        invalid_region2.err().unwrap().kind,
        IndexErrorKind::InvalidRegion
    );

    // fai not exist
    let fai_reader = IndexedReader::from_path("tests/data/fai_no_exist.fa");
    assert_eq!(fai_reader.err().unwrap().kind, IndexErrorKind::FaiIo);

    // fai format error
    let fai_reader = IndexedReader::from_path("tests/data/bad_fai.fa");
    assert_eq!(
        fai_reader.err().unwrap().kind,
        IndexErrorKind::FaiFormatError
    );
}
