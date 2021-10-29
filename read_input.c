//
// Created by ixiaohu on 2021/10/29.
//

#include <stdlib.h>
#include <string.h>

#include "cstl/kstring.h"
#include "read_input.h"


int is_reads_file(const char *fn) {
	gzFile f = gzopen(fn, "r"); assert(f != NULL);
	char buf[2048]; memset(buf, 0, sizeof(buf));
	gzread(f, buf, sizeof(buf));
	gzclose(f);
	int i;
	for (i = 0; i < 2048; i++) {
		if (buf[i] == '@' || buf[i] == '>') return 0; // A FASTQ/FASTA file
	}
	return 1;
}

bseq1_t *load_reads(long n_has_input, int chunk_size, int *n_, gzFile f) {
	int size = 0, m, n;
	char buf[2048];
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (gzgets(f, buf, sizeof(buf)) != NULL) {
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		int len = strlen(buf);
		if (buf[len-1] != '\n') {
			fprintf(stderr, "[%s] Read length exceeds maximum length 2048\n", __func__ );
			abort();
		} else buf[--len] = '\0';
		kstring_t temp; memset(&temp, 0, sizeof(temp));
		ksprintf(&temp, "Reordered_%ld", n_has_input + n + 1);
		seqs[n].id = n;
		seqs[n].name = temp.s;
		seqs[n].comment = NULL;
		seqs[n].seq = strdup(buf);
		seqs[n].l_seq = len;
		seqs[n].qual = NULL;
		seqs[n].sam = NULL;
		size += seqs[n++].l_seq;
		if (size >= chunk_size && (n&1) == 0) break;
	}
	*n_ = n;
	return seqs;
}