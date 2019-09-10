#[macro_use]
extern crate criterion;
extern crate needletail;

use criterion::Criterion;
use needletail::parse_sequence_reader;
use needletail::sequence::Sequence;
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
            parse_sequence_reader(
                fasta_data,
                |_| {},
                |rec| {
                    let seq = rec.seq.normalize(true);
                    let rc = seq.reverse_complement();
                    for (_, _kmer, was_rc) in seq.canonical_kmers(ksize, &rc) {
                        if !was_rc {
                            n_canonical += 1;
                        }
                        n_total += 1;
                    }
                },
            )
            .unwrap();
            assert_eq!(718_007, n_total);
            assert_eq!(350_983, n_canonical);
        });
    });

    group.bench_function("Bitkmer", |bench| {
        bench.iter(|| {
            let mut n_total = 0;
            let mut n_canonical = 0;
            let fasta_data = Cursor::new(data.clone());
            parse_sequence_reader(
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
            assert_eq!(718_007, n_total);
            assert_eq!(350_983, n_canonical);
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
            assert_eq!(250_000, n_bases);
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
            assert_eq!(250_000, n_bases);
        });
    });

    group.bench_function("Needletail", |bench| {
        bench.iter(|| {
            let fastq_data = Cursor::new(data.clone());
            let mut n_bases = 0;
            parse_sequence_reader(
                fastq_data,
                |_| {},
                |seq| {
                    n_bases += seq.seq.len();
                },
            )
            .unwrap();
            assert_eq!(250_000, n_bases);
        });
    });

    group.bench_function("Needletail (No Buffer)", |bench| {
        use needletail::formats::{FastqParser, RecParser};
        bench.iter(|| {
            let mut reader = FastqParser::from_buffer(&data, true);
            let mut n_bases = 0;
            for seq in reader.by_ref() {
                n_bases += seq.unwrap().seq.len();
            }
            assert_eq!(250_000, n_bases);
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
            assert_eq!(738_580, n_bases);
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
            assert_eq!(738_580, n_bases);
        });
    });

    group.bench_function("Needletail", |bench| {
        bench.iter(|| {
            let fasta_data = Cursor::new(data.clone());
            let mut n_bases = 0;
            parse_sequence_reader(
                fasta_data,
                |_| {},
                |seq| {
                    n_bases += seq.seq.len();
                },
            )
            .unwrap();
            assert_eq!(738_580, n_bases);
        });
    });

    group.bench_function("Needletail (No Buffer)", |bench| {
        use needletail::formats::{FastaParser, RecParser};
        bench.iter(|| {
            let mut reader = FastaParser::from_buffer(&data, true);
            let mut n_bases = 0;
            for rec in reader.by_ref() {
                n_bases += rec.unwrap().seq.strip_returns().len();
            }
            assert_eq!(738_580, n_bases);
        });
    });
}

criterion_group!(io, bench_fasta_file, bench_fastq_file);

criterion_main!(kmers, io);
