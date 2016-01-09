#!/usr/bin/env python
import json
import subprocess

fasta_headers = {
    '>0.1 0 0 7644 f': ('assembler', 'ABySS', 'contig'),
    '>0 7644 126659': ('assembler', 'ABySS', 'scaffold'),
    '>2.1 2 0 312 f': ('assembler', 'ABySS2', 'contig'),
    '>2 312 3408': ('assembler', 'ABySS2', 'scaffold'),
    '>contig_0': ('assembler', 'Allpaths-LG', 'contig'),
    '>scaffold_0': ('assembler', 'Allpaths-LG', 'scaffold'),
    '>scf1_0': ('assembler', 'Bambus2', 'contig'),
    '>scf1': ('assembler', 'Bambus2', 'scaffold'),
    '>ctg7180000000867': ('assembler', 'MSR-CA', 'contig'),
    '>scf7180000001028': ('assembler', 'MSR-CA', 'scaffold'),
    '>ctg220002734983': ('assembler', 'CABOG', 'contig'),
    '>scf220002735413': ('assembler', 'CABOG', 'scaffold'),
    '>contig-5594 101 0': ('assembler', 'SGA', 'contig'),
    '>scaffold-0 289': ('assembler', 'SGA', 'scaffold'),
    '>C1742.1 C1742 0 101 f': ('assembler', 'SOAPdenovo', 'contig'),
    '>C1742  1.0': ('assembler', 'SOAPdenovo', 'scaffold'),
    '>velvet.1.1 1 1488': ('assembler', 'Velvet', 'contig'),
    '>velvet.1 NODE_1_length_1458_cov_60.585735': ('assembler', 'Velvet', 'scaffold'),
    '>gi|566658894|gb|AYMX02000001.1| Escherichia coli ATCC BAA-2193 Cont1,whole genome shotgun sequence': ('db', 'NCBI'),  # noqa
    '>ENA|BN000065|BN000065.1 TPA: Homo sapiens SMP1 gene, RHD gene and RHCE gene': ('db', 'European Nucleotide Archive'),  # noqa
    '>gi|54109743|emb|CAG17417.1| Histone [Cotesia congregata bracovirus]': ('db', 'NCBI Protein'),  # from GenBank?  # noqa
    '>gnl|SRA|ERR834220.1 channel_119_read_2': ('db', 'NCBI Sequence Read Archive'),  # from nanopore  # noqa
    '@NS500457:8:H05C5AFXX:1:11101:23206:1045 1:N:0:ATAGCGG+CTATTCC': ('sequencer', 'illumina', ''),
    '@HWUSI-EAS100R:6:73:941:1973#0/1': ('sequencer', 'illumina', ''),
}


def _run_sniffer(executable, check_fasta):
    p = subprocess.Popen([executable], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    print p.communicate(check_fasta)[0]
    return json.loads(p.communicate(check_fasta)[0])


def test_sniffer(executable):
    resp = _run_sniffer(executable, '>test\nAGCT')

    assert resp['compression'] == 'none'
    assert resp['file_type'] == 'fasta'
    assert resp['seq_type'] == 'dna'

    assert not resp['seq_multiline']
    assert not resp['seq_has_gaps']
    assert not resp['seq_has_lowercase']
    assert not resp['seq_has_iupac']
    assert not resp['seq_has_unknowns']

    assert resp['seq_est_avg_len'] == 4
    assert resp['seq_est_gc'] == 0.5

    assert not resp['interleaved']

    resp = _run_sniffer(executable, '@test\nAGCU\n+test\nAAAA')

    assert resp['compression'] == 'none'
    assert resp['file_type'] == 'fastq'
    assert resp['seq_type'] == 'rna'


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Unit tester for FASTX sniffers')
    parser.add_argument('executable', default='./sniffer.py', nargs='?', help='Program to test')

    args = parser.parse_args()
    test_sniffer(args.executable)
