//
// Created by ixiaohu on 2021/9/21.
//

#ifndef ZIP_SEEDING_ZIPMEM_H
#define ZIP_SEEDING_ZIPMEM_H

#include "bwamem/bwamem.h"

#ifdef __cplusplus
extern "C" {
#endif

void zip_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac,
                      int64_t n_processed, int n, bseq1_t *seqs);

#ifdef __cplusplus
}
#endif

#endif //ZIP_SEEDING_ZIPMEM_H
