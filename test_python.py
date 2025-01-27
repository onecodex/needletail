import unittest
from pathlib import Path

from needletail import (
    NeedletailError,
    PyFastxReader,
    Record,
    normalize_seq,
    parse_fastx_file,
    parse_fastx_string,
    reverse_complement,
)

FASTA_FILE = "./tests/data/test.fa"
FASTQ_FILE = "./tests/specimen/FASTQ/example.fastq"


class RecordClassTestCase(unittest.TestCase):
    def test_fasta_record(self):
        record = Record("test description", "AGCTGATCGA")
        self.assertEqual(record.id, "test description")
        self.assertEqual(record.seq, "AGCTGATCGA")
        self.assertIsNone(record.qual)

    def test_fastq_record(self):
        record = Record("test description", "AGCTGATCGA", ";**9;;????")
        self.assertEqual(record.id, "test description")
        self.assertEqual(record.seq, "AGCTGATCGA")
        self.assertEqual(record.qual, ";**9;;????")

    def test_record_properties(self):
        record = Record("test description", "AGCTGATCGA")
        self.assertEqual(record.name, "test")
        self.assertEqual(record.description, "description")

    def test_record_normalize(self):
        record = Record("test", "AGCTGYrtcga")
        record.normalize(iupac=True)
        self.assertEqual(record.seq, "AGCTGYRTCGA")
        record.normalize()
        self.assertEqual(record.seq, "AGCTGNNTCGA")

    def test_format_record_method(self):
        record = Record("test", "AGCTGATCGA")
        self.assertTrue(record.is_fasta())
        self.assertFalse(record.is_fastq())
        record = Record("test", "AGCTGATCGA", ";**9;;????")
        self.assertFalse(record.is_fasta())
        self.assertTrue(record.is_fastq())

    def test_record_eq(self):
        record1 = Record("test", "AGCTGATCGA", ";**9;;????")
        record2 = Record("test", "AGCTGATCGA", ";**9;;????")
        record3 = Record("test2", "AGCTGATCGA", ";**9;;????")
        record4 = Record("test", "TCGATCAGCT", ";**9;;????")
        record5 = Record("test", "AGCTGATCGA", "????;**9;;")
        record6 = Record("test", "AGCTGATCGA")
        self.assertEqual(record1, record2)
        self.assertNotEqual(record1, record3)
        self.assertNotEqual(record1, record4)
        self.assertNotEqual(record1, record5)
        self.assertNotEqual(record1, record6)

    def test_record_str(self):
        self.assertEqual(str(Record("test", "AGCTGATCGA")), ">test\nAGCTGATCGA\n")
        self.assertEqual(
            str(Record("test", "AGCTGATCGA", ";**9;;????")),
            "@test\nAGCTGATCGA\n+\n;**9;;????\n",
        )

    def test_record_repr(self):
        self.assertEqual(
            repr(Record("test", "AGCTGATCGAAGCTGATCGAA")),
            "Record(id=test, seq=AGCTGATCGAAGCTGA…GAA, qual=None)",
        )
        self.assertEqual(
            repr(Record("test", "AGCTGATCGAAGCTGATCGAA", ";**9;;????;**9;;????;")),
            "Record(id=test, seq=AGCTGATCGAAGCTGA…GAA, qual=;**9;;????;**9;;…??;)",
        )

    def test_record_len(self):
        self.assertEqual(len(Record("test", "AGCTGATCGA")), 10)

    def test_record_hash(self):
        record1 = Record("test", "AGCTGATCGA")
        record2 = Record("test", "AGCTGATCGA")
        record3 = Record("test", "AGCTGATCGA", ";**9;;????")
        record4 = Record("test", "AGCTGATCGA", ";**9;;????")
        record5 = Record("test", "TCGATCAGCT")
        record6 = Record("test2", "AGCTGATCGA")
        record7 = Record("test", "AGCTGATCGA", "????;**9;;")
        self.assertEqual(hash(record1), hash(record2))
        self.assertNotEqual(hash(record1), hash(record3))
        self.assertNotEqual(hash(record1), hash(record5))
        self.assertNotEqual(hash(record1), hash(record6))
        self.assertNotEqual(hash(record1), hash(record3))
        self.assertEqual(hash(record3), hash(record4))
        self.assertNotEqual(hash(record3), hash(record7))


class NormalizeTestCase(unittest.TestCase):
    def test_no_normalization_needed(self):
        self.assertEqual(normalize_seq("ACGTU", iupac=False), "ACGTT")

    def test_capitalization(self):
        self.assertEqual(normalize_seq("acgtu", iupac=False), "ACGTT")

    def test_default_parameters(self):
        self.assertEqual(
            normalize_seq("BDHVRYSWKM"), normalize_seq("BDHVRYSWKM", iupac=False)
        )

    def test_iupac_parameter(self):
        self.assertEqual(normalize_seq("BDHVRYSWKM", iupac=False), "NNNNNNNNNN")
        self.assertEqual(normalize_seq("BDHVRYSWKM", iupac=True), "BDHVRYSWKM")
        self.assertEqual(normalize_seq("bdhvryswkm", iupac=True), "BDHVRYSWKM")

    def test_gap_normalization(self):
        self.assertEqual(normalize_seq("N-N-N-N", iupac=False), "N-N-N-N")
        self.assertEqual(normalize_seq("N.N.N.N", iupac=False), "N-N-N-N")
        self.assertEqual(normalize_seq("N~N~N~N", iupac=False), "N-N-N-N")

    def test_whitespace_removal(self):
        self.assertEqual(normalize_seq("N N N N", iupac=False), "NNNN")
        self.assertEqual(normalize_seq("N\tN\tN\tN", iupac=False), "NNNN")
        self.assertEqual(normalize_seq("N\nN\nN\nN", iupac=False), "NNNN")
        self.assertEqual(normalize_seq("N\rN\rN\rN", iupac=False), "NNNN")

    def test_non_alphabet_characters(self):
        self.assertEqual(normalize_seq("N!N!N!N", iupac=False), "NNNNNNN")
        self.assertEqual(normalize_seq("N@N@N@N", iupac=False), "NNNNNNN")
        self.assertEqual(normalize_seq("N#N#N#N", iupac=False), "NNNNNNN")
        self.assertEqual(normalize_seq("N$N$N$N", iupac=False), "NNNNNNN")
        self.assertEqual(normalize_seq("N%N%N%N", iupac=False), "NNNNNNN")
        self.assertEqual(normalize_seq("N^N^N^N", iupac=False), "NNNNNNN")
        self.assertEqual(normalize_seq("N&N&N&N", iupac=False), "NNNNNNN")
        self.assertEqual(normalize_seq("N*N*N*N", iupac=False), "NNNNNNN")
        self.assertEqual(normalize_seq("N|N|N|N", iupac=False), "NNNNNNN")
        self.assertEqual(normalize_seq("N9N5N1N", iupac=False), "NNNNNNN")


class ReverseComplementTestCase(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("a"), "t")
        self.assertEqual(reverse_complement("c"), "g")
        self.assertEqual(reverse_complement("g"), "c")
        self.assertEqual(reverse_complement("n"), "n")
        self.assertEqual(reverse_complement("atcg"), "cgat")


class FileParsingTestCase(unittest.TestCase):
    def get_fasta_reader(self):
        return parse_fastx_file(FASTA_FILE)

    def get_fastq_reader(self):
        return parse_fastx_file(FASTQ_FILE)

    def test_can_parse_fasta_file(self):
        for i, record in enumerate(self.get_fasta_reader()):
            if i == 0:
                self.assertEqual(record.id, "test")
                self.assertEqual(record.seq, "AGCTGATCGA")
                self.assertIsNone(record.qual)
            if i == 1:
                self.assertEqual(record.id, "test2")
                self.assertEqual(record.seq, "TAGC")
                self.assertIsNone(record.qual)
            self.assertTrue(i <= 1)

    def test_can_parse_fastq_file(self):
        for i, record in enumerate(self.get_fastq_reader()):
            if i == 0:
                self.assertEqual(record.id, "EAS54_6_R1_2_1_413_324")
                self.assertEqual(record.seq, "CCCTTCTTGTCTTCAGCGTTTCTCC")
                self.assertEqual(record.qual, ";;3;;;;;;;;;;;;7;;;;;;;88")
            if i == 1:
                self.assertEqual(record.id, "EAS54_6_R1_2_1_540_792")
                self.assertEqual(record.seq, "TTGGCAGGCCAAGGCCGATGGATCA")
                self.assertEqual(record.qual, ";;;;;;;;;;;7;;;;;-;;;3;83")
            self.assertTrue(i <= 2)

    def test_pathlib_path_input(self):
        self.assertIsInstance(parse_fastx_file(Path(FASTA_FILE)), PyFastxReader)


class StrParsingTestCase(FileParsingTestCase):
    def get_fasta_reader(self):
        with open(FASTA_FILE) as f:
            content = f.read()
            return parse_fastx_string(content)

    def get_fastq_reader(self):
        with open(FASTQ_FILE) as f:
            content = f.read()
            return parse_fastx_string(content)

    def test_pathlib_path_input(self):
        pass


class ErroringTestCase(unittest.TestCase):
    def test_file_not_found(self):
        with self.assertRaises(NeedletailError):
            parse_fastx_file("hey")

    def test_invalid_record(self):
        with self.assertRaises(NeedletailError):
            for i, record in enumerate(parse_fastx_string("Not a valid file")):
                print(i)


if __name__ == "__main__":
    unittest.main()
