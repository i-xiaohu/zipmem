//
// Created by ixiaohu on 2021/9/21.
//

#ifndef ZIP_SEEDING_ZIPMEM_H
#define ZIP_SEEDING_ZIPMEM_H

#include "bwamem/bwamem.h"

typedef struct {
	int n; // consensus sequence length
	uint8_t *s; // consensus sequence
	int *cnt[4]; // Saving votes.
} con_seq_t;

typedef struct {
	size_t n, m;
	con_seq_t *a;
} con_seq_v;

// Allocated memory for all CS.
typedef struct {
	uint8_t *s; // CS bases
	int *cnt[4]; // counts of A,C,G,T in CS
	int *mem_loc; // memory location for initializing CS
	int *offset; // offset on CS of all reads
	uint8_t *is_rc; // if reads are reversed/complement in CS
	int *belong_to; // which CS the read belongs to
} cs_aux_t;

#ifdef __cplusplus
extern "C" {
#endif
	cs_aux_t* cs_aux_init(int n, const bseq1_t *reads);
	con_seq_v connect_reads_to_cs(int threads_n, int n, bseq1_t *reads, cs_aux_t *aux);
	void cs_aux_destroy(cs_aux_t *a);

	/**
	 * Seeding for reordered reads
	 * Here, $seqs[i].sam is storing seeds in format of seed clusters.
	 * Specifically, (qb, l, score, occ, pos1, pos2, pos3 ...)
	 * All seeds are encoded and outputted to a binary file.
	 * @param opt Seeding options: seed length, minimum occurrence, etc.
	 * @param bwt FM-index
	 * @param n Reads number
	 * @param seqs Reordered reads
	 */
	void zipmem_seeding(const mem_opt_t *opt, const bwt_t *bwt, int n, bseq1_t *seqs);

#ifdef __cplusplus
}
#endif

#endif //ZIP_SEEDING_ZIPMEM_H
