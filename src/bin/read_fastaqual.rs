use needletail::parse_fastaqual_file;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let fasta_path = &args[1];
    let qual_path = &args[2];
    let mut reader = parse_fastaqual_file(fasta_path, qual_path).expect("invalid reader");
    let mut record_num = 1;

    while let Some(r) = reader.next() {
        let record = r.expect("invalid record");
        let id = std::str::from_utf8(record.id()).expect("invalid id");
        let num_bases = record.num_bases();
        let scores = record.qual().expect("invalid quality scores");
        let min_score = scores.iter().min().expect("record is empty");
        let max_score = scores.iter().max().expect("record is empty");

        println!(
            "Record {} has id '{}', {} bases, and quality scores are in the range [{}, {}]",
            record_num, id, num_bases, min_score, max_score
        );
        record_num += 1;
    }
}
