//
// Created by ixiaohu on 2021/10/29.
//

#ifndef ZIP_SEEDING_READ_INPUT_H
#define ZIP_SEEDING_READ_INPUT_H

#include "bwalib/bwa.h"

#ifdef __cplusplus
extern "C" {
#endif
	int is_reads_file(const char *fn);
	bseq1_t *load_reads(long n_has_input, int chunk_size, int *n_, gzFile f);
#ifdef __cplusplus
}
#endif


#endif //ZIP_SEEDING_READ_INPUT_H
