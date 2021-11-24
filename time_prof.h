//
// Created by ixiaohu on 2020/5/31.
//

#ifndef ZIP_SEEDING_5_9_TIME_PROF_H
#define ZIP_SEEDING_5_9_TIME_PROF_H

typedef struct {
	long reads_n;
	long read_len;
	long cs_bases;
	long mis_n[256];
	long reused_n[256];
	long sup_n[256];
	// Time
	double t_gencs[2];
	double t_csseed[2];
	double t_seeding[2];
} zsmem_prof_t;

typedef struct {
	double t_seeding[2];
	double t_extending[2];
} bwamem_prof_t;

static void zsmem_prof_output(zsmem_prof_t *a) {
	int i;
	for (i = 1; i < 256; i++) a->mis_n[0] += a->mis_n[i];
	for (i = 1; i < 256; i++) a->reused_n[0] += a->reused_n[i];
	for (i = 1; i < 256; i++) a->sup_n[0] += a->sup_n[i];
	fprintf(stderr, "Reads number        \t%ld\n", a->reads_n);
	fprintf(stderr, "Read length         \t%ld\n", a->read_len);
	fprintf(stderr, "Compression ratio   \t%.2f X\n", 1.0*(a->reads_n*a->read_len)/a->cs_bases);
	fprintf(stderr, "Mismatch percentage \t%.2f %%\n", 100.0*a->mis_n[0]/(a->reads_n*a->read_len));
	fprintf(stderr, "Reused seeds        \t%.2f\n", 1.0*a->reused_n[0]/a->reads_n);
	fprintf(stderr, "Supplementary seeds \t%.2f\n", 1.0*a->sup_n[0]/a->reads_n);
}

#endif //ZIP_SEEDING_5_9_TIME_PROF_H
