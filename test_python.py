import unittest

from needletail import parse_fastx_file, NeedletailError, reverse_complement, normalize_seq


class ParsingTestCase(unittest.TestCase):
    def test_can_parse_fasta(self):
        for i, record in enumerate(parse_fastx_file("./tests/data/test.fa")):
            if i == 0:
                self.assertEqual(record.id, "test")
                self.assertEqual(record.seq, "AGCTGATCGA")
                self.assertIsNone(record.qual)
                record.normalize(iupac=False)
                self.assertEqual(record.seq, "AGCTGATCGA")
                self.assertTrue(record.is_fasta())
            if i == 1:
                self.assertEqual(record.id, "test2")
                self.assertEqual(record.seq, "TAGC")
                self.assertIsNone(record.qual)
                record.normalize(iupac=False)
                self.assertEqual(record.seq, "TAGC")
                self.assertTrue(record.is_fasta())

            self.assertTrue(i <= 1)

    def test_can_parse_fastq(self):
        for i, record in enumerate(parse_fastx_file("./tests/specimen/FASTQ/example.fastq")):
            if i == 0:
                self.assertEqual(record.id, "EAS54_6_R1_2_1_413_324")
                self.assertEqual(record.seq, "CCCTTCTTGTCTTCAGCGTTTCTCC")
                self.assertEqual(record.qual, ";;3;;;;;;;;;;;;7;;;;;;;88")
                record.normalize(iupac=False)
                self.assertEqual(record.seq, "CCCTTCTTGTCTTCAGCGTTTCTCC")
                self.assertTrue(record.is_fastq())
            if i == 1:
                self.assertEqual(record.id, "EAS54_6_R1_2_1_540_792")
                self.assertEqual(record.seq, "TTGGCAGGCCAAGGCCGATGGATCA")
                self.assertEqual(record.qual, ";;;;;;;;;;;7;;;;;-;;;3;83")
                record.normalize(iupac=False)
                self.assertEqual(record.seq, "TTGGCAGGCCAAGGCCGATGGATCA")
                self.assertTrue(record.is_fastq())

            self.assertTrue(i <= 2)


class MiscelleanousTestCase(unittest.TestCase):
    def test_normalize_seq(self):
        self.assertEqual(normalize_seq("ACGTU", iupac=False), "ACGTT")
        self.assertEqual(normalize_seq("acgtu", iupac=False), "ACGTT")
        self.assertEqual(normalize_seq("N.N-N~N N", iupac=False), "N-N-N-NN")
        self.assertEqual(normalize_seq("BDHVRYSWKM", iupac=True), "BDHVRYSWKM")
        self.assertEqual(normalize_seq("bdhvryswkm", iupac=True), "BDHVRYSWKM")

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("a"), "t")
        self.assertEqual(reverse_complement("c"), "g")
        self.assertEqual(reverse_complement("g"), "c")
        self.assertEqual(reverse_complement("n"), "n")

        self.assertEqual(reverse_complement("atcg"), "cgat")


if __name__ == '__main__':
    unittest.main()
