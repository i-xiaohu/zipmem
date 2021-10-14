//
// Created by ixiaohu on 2021/10/14.
//

#include <stdio.h>

#include "bwalib/bwa.h"
#include "bwalib/kseq.h"

KSEQ_DECLARE(gzFile)

void check_reused_seeds(const char *idx_name, const char *fastq_fn, const char *seed_fn) {
	bwaidx_t *idx = bwa_idx_load_from_disk(idx_name, BWA_IDX_ALL);
	gzFile fp = gzopen(fastq_fn, "r");
	kseq_t *ks = kseq_init(fp);
	FILE *fseeds = fopen(seed_fn, "r");
	long reads_n = 0;
	long seeds_n = 0, no_mismatch = 0, mismatches_n = 0;
	while (1) {
		int i, j, k, h, n_seqs;
		long size = 0;
		bseq1_t *seqs = bseq_read(240 * 1000 * 1000, &n_seqs, ks, NULL);
		if (seqs == NULL) break;
		reads_n += n_seqs;
		for (i = 0; i < n_seqs; ++i) size += seqs[i].l_seq;
		long total_bytes = 0;
		for (i = 0; i < n_seqs; i++) {
			bseq1_t *read = &seqs[i];
			for (j = 0; j < read->l_seq; j++) read->seq[j] = nst_nt4_table[read->seq[j]];
			long seed_bytes = 0;
			fread(&seed_bytes, sizeof(long), 1, fseeds);
			total_bytes += seed_bytes;
			uint8_t *seeds = malloc(seed_bytes);
			fread(seeds, sizeof(uint8_t), seed_bytes, fseeds);
			long pointer = 0;
			int mem_n = *(int*)(seeds + pointer); pointer += sizeof(int);
			for (j = 0; j < mem_n; j++) {
				int sb = *(int*)(seeds + pointer); pointer += sizeof(int);
				int sl = *(int*)(seeds + pointer); pointer += sizeof(int);
				long sa_size = *(long*)(seeds + pointer); pointer += sizeof(long);
				int true_occ = sa_size < 500 ?sa_size :500;
				seeds_n += true_occ;
				for (k = 0; k < true_occ; k++) {
					long rb = *(long*)(seeds + pointer); pointer += sizeof(long);
					int rid = bns_intv2rid(idx->bns, rb, rb + sl);
					if (rid < 0) continue;
					int64_t get_len = 0;
					uint8_t *ref = bns_get_seq(idx->bns->l_pac, idx->pac, rb, rb + sl, &get_len);
					if (get_len != sl) {
						fprintf(stderr, "%d %ld\n", sl, get_len);
					}
					assert(get_len == sl);
					int mis = 0;
					for (h = 0; h < sl; h++) {
						if (read->seq[sb + h] != ref[h]) mis++;
					}
					if (!mis) no_mismatch++;
					else mismatches_n += mis;
				}
			}
			assert(pointer == seed_bytes);
			free(seeds);
			free(read->name); free(read->comment); free(read->seq); free(read->qual);
		}
		fprintf(stderr, "[%s] Processed %d reads (%ld Mbp) and %ld MB seeds\n",
		        __func__, n_seqs, (long)size/1000/1000, total_bytes/1024/1024);
	}
	fprintf(stderr, "reads number\t%ld\n", reads_n);
	fprintf(stderr, "seeds number\t%.2f\n", 1.0 * seeds_n / reads_n);
	fprintf(stderr, "Exactly matched seeds\t%.2f\n", 100.0 * no_mismatch / seeds_n);
	fprintf(stderr, "Mismatched number of noisy seeds\t%.2f\n", 1.0 * mismatches_n / (seeds_n-no_mismatch));
}

int main(int argc, char *argv[]) {
	check_reused_seeds(
			"/vol1/agis/ruanjue_group/jifahu/ecoli/ref/ecoli.fa",
			"/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/spring/SRR1562082_1.fq.gz",
			"/vol1/agis/ruanjue_group/jifahu/zsmem-experiments/zipmem/spring_SRR1562082_1.seeds"
	);
}
