#[macro_use]
extern crate bencher;
extern crate needletail;

use bencher::Bencher;
use needletail::fastx;
use std::fs::File;
use std::io::{Cursor, Read};

// from Bio.SeqIO import parse
// n_total = sum([len([k for k in slid_win(i.seq, 31) if set(k).issubset({'A', 'C', 'G', 'T'})]) for i in SeqIO.parse('./tests/data/28S.fasta', 'fasta')])

fn bench_kmer_speed(bench: &mut Bencher) {
    let filename = "tests/data/28S.fasta";
    let ksize = 31;

    bench.iter(|| {
        let mut n_total = 0;
        let mut n_canonical = 0;
        let file = File::open(filename).unwrap();
        fastx::parse_sequences(
            file,
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
    })
}

fn bench_bitkmer_speed(bench: &mut Bencher) {
    let filename = "tests/data/28S.fasta";
    let ksize = 31;

    bench.iter(|| {
        let mut n_total = 0;
        let mut n_canonical = 0;
        let file = File::open(filename).unwrap();
        fastx::parse_sequences(
            file,
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
    })
}

benchmark_group!(kmers, bench_kmer_speed, bench_bitkmer_speed);

fn bench_fastq_bytes(bench: &mut Bencher) {
    let mut data: Vec<u8> = vec![];
    let mut f = File::open("tests/data/PRJNA271013_head.fq").unwrap();
    let _ = f.read_to_end(&mut data);

    bench.iter(|| {
        let mut n_bases = 0;
        fastx::parse_sequences(
            Cursor::new(&data),
            |_| {},
            |seq| {
                n_bases += seq.seq.len();
            },
        )
        .unwrap();
        assert_eq!(250000, n_bases);
    })
}

fn bench_fastq_file(bench: &mut Bencher) {
    let filename = "tests/data/PRJNA271013_head.fq";

    // warming up the cache doesn't seem to make the timings more repeatable?
    // fastx::fastx_file(&filename[..], |seq| { assert!(seq.1.len() > 0) }).unwrap();
    bench.iter(|| {
        let mut n_bases = 0;
        fastx::parse_sequences(
            File::open(filename).unwrap(),
            |_| {},
            |seq| {
                n_bases += seq.seq.len();
            },
        )
        .unwrap();
        assert_eq!(250000, n_bases);
    })
}

fn bench_fasta_bytes(bench: &mut Bencher) {
    let filename = String::from("tests/data/28S.fasta");

    let mut data: Vec<u8> = vec![];
    let mut f = File::open(filename).unwrap();
    let _ = f.read_to_end(&mut data);

    bench.iter(|| {
        let mut n_bases = 0;
        fastx::parse_sequences(
            Cursor::new(&data),
            |_| {},
            |seq| {
                n_bases += seq.seq.len();
            },
        )
        .unwrap();
        assert_eq!(738580, n_bases);
    })
}

fn bench_fasta_file(bench: &mut Bencher) {
    let filename = "tests/data/28S.fasta";

    // warming up the cache doesn't seem to make the timings more repeatable?
    // fastx::fastx_file(&filename[..], |seq| { assert!(seq.1.len() > 0) }).unwrap();
    bench.iter(|| {
        let mut n_bases = 0;
        fastx::parse_sequences(
            File::open(filename).unwrap(),
            |_| {},
            |seq| {
                n_bases += seq.seq.len();
            },
        )
        .unwrap();
        assert_eq!(738580, n_bases);
    })
}

benchmark_group!(
    fastx,
    bench_fastq_bytes,
    bench_fastq_file,
    bench_fasta_bytes,
    bench_fasta_file
);
benchmark_main!(kmers, fastx);
