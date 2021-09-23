//
// Created by ixiaohu on 2020/5/31.
//

#ifndef ZIP_SEEDING_5_9_TIME_PROF_H
#define ZIP_SEEDING_5_9_TIME_PROF_H

#define kind_reused_mem 0
#define kind_reused_rem 1
#define kind_reused_lem 2
#define kind_sup_mem    3
#define kind_sup_rem    4
#define kind_sup_lem    5
#define kind_sig_mem    6
#define kind_sig_rem    7
#define kind_sig_lem    8

typedef struct {
	long reads_n;
	long read_len;
	long singletons;
	long cs_bases;
	// Seeds
	long seeds, seeds_mt[256];
	// Reused MEM, Reused REM, reused LEM, sup MEM, sup REM, sup LEM, sig MEM, sig REM, sig LEM.
	long sub_seeds[9], sub_seeds_mt[9][256];
	long correct_seeds, correct_seeds_mt[256];
	long correct_sub_seeds[9], correct_sub_seeds_mt[9][256];
	long chains, chains_mt[256];
	long regs, regs_mt[256];
	long patch_regs, patch_regs_mt[256];
	long reg2aln, reg2aln_mt[256];
	// Time
	double t_gencs[2];
	double t_csseed[2];
	double t_worker1[2];
	double t_worker2[2];
} zsmem_prof_t;

static void zsmem_prof_sumup(int n, zsmem_prof_t *a) {
	int i, j;
	for (i = 0; i < n; i++) {
		a->seeds += a->seeds_mt[i];
		a->correct_seeds += a->correct_seeds_mt[i];
		for(j = 0; j < 9; j++) {
			a->sub_seeds[j] += a->sub_seeds_mt[j][i];
			a->correct_sub_seeds[j] += a->correct_sub_seeds_mt[j][i];
		}
		a->chains += a->chains_mt[i];
		a->regs += a->regs_mt[i];
		a->patch_regs += a->patch_regs_mt[i];
		a->reg2aln += a->reg2aln_mt[i];
	}
}

static void zsmem_prof_output(const zsmem_prof_t *a) {
	fprintf(stderr, "gencs_ctime %.2f\n", a->t_gencs[0]);
	fprintf(stderr, "gencs_rtime %.2f\n", a->t_gencs[1]);
	fprintf(stderr, "csseed_ctime %.2f\n", a->t_csseed[0]);
	fprintf(stderr, "csseed_rtime %.2f\n", a->t_csseed[1]);
	fprintf(stderr, "worker1_ctime %.2f\n", a->t_worker1[0]);
	fprintf(stderr, "worker1_rtime %.2f\n", a->t_worker1[1]);
	fprintf(stderr, "worker2_ctime %.2f\n", a->t_worker2[0]);
	fprintf(stderr, "worker2_rtime %.2f\n", a->t_worker2[1]);

	fprintf(stderr, "input_reads %ld\n", a->reads_n);
	fprintf(stderr, "read_len %ld\n", a->read_len);
	fprintf(stderr, "singletons %ld\n", a->singletons);
	fprintf(stderr, "cs_bases %ld\n", a->cs_bases);

	fprintf(stderr, "seeds %ld\n", a->seeds);
	fprintf(stderr, "reused_mem_seeds %ld\n", a->sub_seeds[0]);
	fprintf(stderr, "reused_rem_seeds %ld\n", a->sub_seeds[1]);
	fprintf(stderr, "reused_lem_seeds %ld\n", a->sub_seeds[2]);
	fprintf(stderr, "sup_mem_seeds %ld\n", a->sub_seeds[3]);
	fprintf(stderr, "sup_rem_seeds %ld\n", a->sub_seeds[4]);
	fprintf(stderr, "sup_lem_seeds %ld\n", a->sub_seeds[5]);
	fprintf(stderr, "sig_mem_seeds %ld\n", a->sub_seeds[6]);
	fprintf(stderr, "sig_rem_seeds %ld\n", a->sub_seeds[7]);
	fprintf(stderr, "sig_lem_seeds %ld\n", a->sub_seeds[8]);

	fprintf(stderr, "correct_seeds %ld\n", a->correct_seeds);
	fprintf(stderr, "correct_reused_mem_seeds %ld\n", a->correct_sub_seeds[0]);
	fprintf(stderr, "correct_reused_rem_seeds %ld\n", a->correct_sub_seeds[1]);
	fprintf(stderr, "correct_reused_lem_seeds %ld\n", a->correct_sub_seeds[2]);
	fprintf(stderr, "correct_sup_mem_seeds %ld\n", a->correct_sub_seeds[3]);
	fprintf(stderr, "correct_sup_rem_seeds %ld\n", a->correct_sub_seeds[4]);
	fprintf(stderr, "correct_sup_lem_seeds %ld\n", a->correct_sub_seeds[5]);
	fprintf(stderr, "correct_sig_mem_seeds %ld\n", a->correct_sub_seeds[6]);
	fprintf(stderr, "correct_sig_rem_seeds %ld\n", a->correct_sub_seeds[7]);
	fprintf(stderr, "correct_sig_lem_seeds %ld\n", a->correct_sub_seeds[8]);

	fprintf(stderr, "chains %ld\n", a->chains);
	fprintf(stderr, "regs %ld\n", a->regs);
	fprintf(stderr, "patch_reg %ld\n", a->patch_regs);
	fprintf(stderr, "reg2aln %ld\n", a->reg2aln);
}

#endif //ZIP_SEEDING_5_9_TIME_PROF_H
