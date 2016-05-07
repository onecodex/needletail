// using https://github.com/lh3/readfq for parsing
// derived in part from http://stackoverflow.com/questions/19390245/how-to-parse-a-fasta-file-using-kseq-h
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "./readfq/kseq.h"
KSEQ_INIT(FILE*, read)

int main(int argc, char* argv[]) {
	FILE* fp = fopen(argv[1], "r");
	if (fp == 0) {
		perror("fopen");
		exit(1);
	}
	kseq_t *seq = kseq_init(fileno(fp));

	unsigned long slen = 0;
	while (kseq_read(seq) >= 0) {
		slen += seq->seq.l;
	}

	printf("%lu\n", slen);
	
	kseq_destroy(seq);
	fclose(fp);
	return (0);
}
