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

	fprintf(stderr, "gencs_time  \t%.1f\t%.1f\t%.1fX\n",
		 a->t_gencs[0], a->t_gencs[1], a->t_gencs[0]/a->t_gencs[1]);
	fprintf(stderr, "csseed_time \t%.1f\t%.1f\t%.1fX\n",
		 a->t_csseed[0], a->t_csseed[1], a->t_csseed[0]/a->t_csseed[1]);
	fprintf(stderr, "input_reads\t%ld\n", a->reads_n);
	fprintf(stderr, "read_len\t%ld\n", a->read_len);
	fprintf(stderr, "Compression ration  \t%.2f\n", 1.0*a->reads_n*a->read_len/a->cs_bases);
	fprintf(stderr, "Mismatches number   \t%.2f\n", 1.0*a->mis_n[0]/a->reads_n);
	fprintf(stderr, "Reused seeds        \t%.2f\n", 1.0*a->reused_n[0]/a->reads_n);
	fprintf(stderr, "Supplementary seeds \t%.2f\n", 1.0*a->sup_n[0]/a->reads_n);
}

#endif //ZIP_SEEDING_5_9_TIME_PROF_H
