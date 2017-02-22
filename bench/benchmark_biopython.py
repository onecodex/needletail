from Bio import SeqIO
import sys

total = 0
with open(sys.argv[1], "rU") as handle:
    for record in SeqIO.parse(handle, "fastq"):
        total += len(record.seq)

print(total)
