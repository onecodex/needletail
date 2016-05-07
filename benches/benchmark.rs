#[macro_use]
extern crate bencher;
extern crate needletail;

use needletail::{fastx, kmer, bitkmer};
use bencher::Bencher;

fn check_kmer_speed(bench: &mut Bencher) {
    let filename = String::from("tests/data/28S.fasta");
    let ksize = 31;

    bench.iter(|| {
        let mut n_total = 0;
        let mut n_canonical = 0;
        fastx::fastx_file(&filename[..], |seq| {
            // let normalized_seq = kmer::normalize(Cow::Borrowed(&seq.1), true);
            for k in seq.1.windows(ksize) {
                let l = kmer::canonical(k).into_owned();
                if l == k {
                    n_canonical += 1;
                }
                n_total += 1;
            }
        });
        assert_eq!(213703, n_total);
        assert_eq!(108521, n_canonical);
    })
}

fn check_bitkmer_speed(bench: &mut Bencher) {
    let filename = String::from("tests/data/28S.fasta");
    let ksize = 31;

    bench.iter(|| {
        let mut n_total = 0;
        let mut n_canonical = 0;
        fastx::fastx_file(&filename[..], |seq| {
            for k in bitkmer::BitKmerIter::new(&seq.1, ksize) {
                let l = bitkmer::canonical(k);
                if l == k {
                    n_canonical += 1;
                }
                n_total += 1;
            }
        });
        assert_eq!(213703, n_total);
        assert_eq!(108521, n_canonical);
    })
}


benchmark_group!(benches, check_kmer_speed, check_bitkmer_speed);
benchmark_main!(benches);
