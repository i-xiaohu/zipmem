#include <stdio.h>
#include "bwalib/bwa.h"
#include "bwalib/kseq.h"
#include "bwalib/kopen.h"
#include "bwalib/utils.h"
#include "bwamem/bwamem.h"
#include "time_prof.h"

KSEQ_DECLARE(gzFile)

zsmem_prof_t zmp;

int main(int argc, char *argv[]) {
	bwaidx_t *idx = bwa_idx_load_from_shm(argv[1]);
	if (idx == NULL) {
		if ((idx = bwa_idx_load(argv[1], BWA_IDX_ALL)) == NULL) return 1;
		else fprintf(stderr, "[%s] Load the FM-index from disk\n", __func__ );
	} else {
		fprintf(stderr, "[%s] Load the FM-index from shared memory\n", __func__);
	}
	int fd;
	void *ko = kopen(argv[2], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[%s] Fail to open file `%s'.\n", __func__, argv[2]);
		return 1;
	}
	gzFile fp = gzdopen(fd, "r");
	kseq_t *ks = kseq_init(fp);
	int actual_chunk_size = 10 * 1000 * 1000;
	int n_seqs;
	mem_opt_t *opt = mem_opt_init();
	bseq1_t *seqs = bseq_read(actual_chunk_size, &n_seqs, ks, NULL);

	mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, 0, n_seqs, seqs, NULL);
	int i;
	for (i = 0; i < n_seqs; i++) {
		fprintf(stdout, "%s", seqs[i].sam);
	}
	bwa_idx_destroy(idx);
	kseq_destroy(ks);
	err_gzclose(fp);
	kclose(ko);
	return 0;
}
