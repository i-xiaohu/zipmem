//
// Created by ixiaohu on 2021/9/21.
//

#ifndef ZIP_SEEDING_ZIPMEM_H
#define ZIP_SEEDING_ZIPMEM_H

#include "bwamem/bwamem.h"

typedef struct {
	int l, r, n; // consensus sequence length
	uint8_t *s; // consensus sequence
	uint8_t *seeds; // CS seeds
} con_seq_t;

typedef struct {
	size_t n, m;
	con_seq_t *a;
} con_seq_v;


#ifdef __cplusplus
extern "C" {
#endif
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
