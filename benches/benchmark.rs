#[macro_use]
extern crate criterion;
extern crate needletail;

use criterion::Criterion;
use needletail::parse_sequences;
use std::fs::File;
use std::io::{Cursor, Read};

// from Bio.SeqIO import parse
// n_total = sum([len([k for k in slid_win(i.seq, 31) if set(k).issubset({'A', 'C', 'G', 'T'})]) for i in SeqIO.parse('./tests/data/28S.fasta', 'fasta')])

fn bench_kmer_speed(c: &mut Criterion) {
    let ksize = 31;

    let mut data: Vec<u8> = vec![];
    let mut f = File::open("tests/data/28S.fasta").unwrap();
    let _ = f.read_to_end(&mut data);

    let mut group = c.benchmark_group("Kmerizing");
    group.sample_size(10);

    group.bench_function("Kmer", |b| {
        b.iter(|| {
            let mut n_total = 0;
            let mut n_canonical = 0;
            let fasta_data = Cursor::new(data.clone());
            parse_sequences(
                fasta_data,
                |_| {},
                |seq| {
                    for (_, _kmer, was_rc) in seq.normalize(true).kmers(ksize, true) {
                        if !was_rc {
                            n_canonical += 1;
                        }
                        n_total += 1;
                    }
                },
            )
            .unwrap();
            assert_eq!(718007, n_total);
            assert_eq!(350983, n_canonical);
        });
    });

    group.bench_function("Bitkmer", |bench| {
        bench.iter(|| {
            let mut n_total = 0;
            let mut n_canonical = 0;
            let fasta_data = Cursor::new(data.clone());
            parse_sequences(
                fasta_data,
                |_| {},
                |seq| {
                    for (_, _kmer, was_rc) in seq.bit_kmers(ksize, true) {
                        if !was_rc {
                            n_canonical += 1;
                        }
                        n_total += 1;
                    }
                },
            )
            .unwrap();
            assert_eq!(718007, n_total);
            assert_eq!(350983, n_canonical);
        });
    });
}

criterion_group!(kmers, bench_kmer_speed);

fn bench_fastq_file(c: &mut Criterion) {
    use bio::io::fastq as bio_fastq;
    use seq_io::fastq as seq_fastq;
    use seq_io::fastq::Record;

    let mut data: Vec<u8> = vec![];
    let mut f = File::open("tests/data/PRJNA271013_head.fq").unwrap();
    let _ = f.read_to_end(&mut data);

    let mut group = c.benchmark_group("FASTQ parsing");

    group.bench_function("RustBio", |bench| {
        bench.iter(|| {
            let fastq_data = Cursor::new(data.clone());
            let reader = bio_fastq::Reader::new(fastq_data);
            let mut n_bases = 0;
            for record in reader.records() {
                let record = record.unwrap();
                n_bases += record.seq().len()
            }
            assert_eq!(250000, n_bases);
        });
    });

    group.bench_function("SeqIO", |bench| {
        bench.iter(|| {
            let fastq_data = Cursor::new(data.clone());
            let mut reader = seq_fastq::Reader::new(fastq_data);
            let mut n_bases = 0;
            while let Some(result) = reader.next() {
                let record = result.unwrap();
                let seqlen = record.seq().len();
                n_bases += seqlen;
            }
            assert_eq!(250000, n_bases);
        });
    });

    group.bench_function("Needletail", |bench| {
        bench.iter(|| {
            let fastq_data = Cursor::new(data.clone());
            let mut n_bases = 0;
            parse_sequences(
                fastq_data,
                |_| {},
                |seq| {
                    n_bases += seq.seq.len();
                },
            )
            .unwrap();
            assert_eq!(250000, n_bases);
        });
    });

    group.bench_function("Needletail (Macro)", |bench| {
        use needletail::formats::FastqReader;
        use needletail::{parse_stream, ParseError};
        #[inline]
        fn get_n_bases(mut fastq_data: &mut dyn Read) -> Result<usize, ParseError> {
            let mut n_bases = 0;
            parse_stream!(&mut fastq_data, &b""[..], FastqReader, rec, {
                n_bases += rec.seq.len();
            });
            Ok(n_bases)
        }

        bench.iter(|| {
            let mut fastq_data = Cursor::new(data.clone());
            let n_bases = get_n_bases(&mut fastq_data).unwrap();
            assert_eq!(250000, n_bases);
        });
    });

    group.bench_function("Needletail (No Buffer)", |bench| {
        use needletail::formats::{FastqReader, RecReader};
        bench.iter(|| {
            let mut reader = FastqReader::from_buffer(&data, true);
            let mut n_bases = 0;
            for seq in reader.by_ref() {
                n_bases += seq.unwrap().seq.len();
            }
            assert_eq!(250000, n_bases);
        });
    });
}

fn bench_fasta_file(c: &mut Criterion) {
    use bio::io::fasta as bio_fasta;
    use seq_io::fasta as seq_fasta;

    let mut data: Vec<u8> = vec![];
    let mut f = File::open("tests/data/28S.fasta").unwrap();
    let _ = f.read_to_end(&mut data);

    let mut group = c.benchmark_group("FASTA parsing");

    group.bench_function("RustBio", |bench| {
        bench.iter(|| {
            let fasta_data = Cursor::new(data.clone());
            let reader = bio_fasta::Reader::new(fasta_data);
            let mut n_bases = 0;
            for record in reader.records() {
                let record = record.unwrap();
                n_bases += record.seq().len()
            }
            assert_eq!(738580, n_bases);
        });
    });

    group.bench_function("SeqIO", |bench| {
        bench.iter(|| {
            let fasta_data = Cursor::new(data.clone());
            let mut reader = seq_fasta::Reader::new(fasta_data);
            let mut n_bases = 0;
            while let Some(result) = reader.next() {
                let record = result.unwrap();
                for s in record.seq_lines() {
                    n_bases += s.len();
                }
            }
            assert_eq!(738580, n_bases);
        });
    });

    group.bench_function("Needletail", |bench| {
        bench.iter(|| {
            let fasta_data = Cursor::new(data.clone());
            let mut n_bases = 0;
            parse_sequences(
                fasta_data,
                |_| {},
                |seq| {
                    n_bases += seq.seq.len();
                },
            )
            .unwrap();
            assert_eq!(738580, n_bases);
        });
    });

    group.bench_function("Needletail (Macro)", |bench| {
        use needletail::formats::FastaReader;
        use needletail::seq::Sequence;
        use needletail::{parse_stream, ParseError};
        #[inline]
        fn get_n_bases(mut fasta_data: &mut dyn Read) -> Result<usize, ParseError> {
            let mut n_bases = 0;
            parse_stream!(&mut fasta_data, &b""[..], FastaReader, rec, {
                let seq = Sequence::from(rec);
                n_bases += seq.seq.len();
            });
            Ok(n_bases)
        }

        bench.iter(|| {
            let mut fasta_data = Cursor::new(data.clone());
            let n_bases = get_n_bases(&mut fasta_data).unwrap();
            assert_eq!(738580, n_bases);
        });
    });

    group.bench_function("Needletail (No Buffer)", |bench| {
        use needletail::formats::{FastaReader, RecReader};
        use needletail::seq::Sequence;
        bench.iter(|| {
            let mut reader = FastaReader::from_buffer(&data, true);
            let mut n_bases = 0;
            for rec in reader.by_ref() {
                let seq = Sequence::from(rec.unwrap());
                n_bases += seq.seq.len();
            }
            assert_eq!(738580, n_bases);
        });
    });
}

criterion_group!(io, bench_fasta_file, bench_fastq_file);

criterion_main!(kmers, io);
