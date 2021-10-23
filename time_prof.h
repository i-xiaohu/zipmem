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
	double input[2];
	double process[2];
	double output[2];
	double allocate[2];
	double assemble[2];
	double sumup[2];
	double t_gencs[2];
	double t_csseed[2];
	double t_supseed[2];
} zsmem_prof_t;

static void zsmem_prof_output(zsmem_prof_t *a) {
	int i;
	for (i = 1; i < 256; i++) a->mis_n[0] += a->mis_n[i];
	for (i = 1; i < 256; i++) a->reused_n[0] += a->reused_n[i];
	for (i = 1; i < 256; i++) a->sup_n[0] += a->sup_n[i];

	fprintf(stderr, "Input_time  \t%.1f\n", a->input[1]);
	fprintf(stderr, "Process_time\t%.1f\t%.1f\t%.1fX\n",
		 a->process[0], a->process[1], a->process[0]/a->process[1]);
	fprintf(stderr, "output_time \t%.1f\n", a->output[1]);
	fprintf(stderr, "gencs_time  \t%.1f\t%.1f\t%.1fX\n",
		 a->t_gencs[0], a->t_gencs[1], a->t_gencs[0]/a->t_gencs[1]);
	fprintf(stderr, "  Allocating\t%.1f\t%.1f\n", a->allocate[0], a->allocate[1]);
	fprintf(stderr, "  Assemble  \t%.1f\t%.1f\n", a->assemble[0], a->assemble[1]);
	fprintf(stderr, "  Sumup     \t%.1f\t%.1f\n", a->sumup[0], a->sumup[1]);
	fprintf(stderr, "csseed_time \t%.1f\t%.1f\t%.1fX\n",
		 a->t_csseed[0], a->t_csseed[1], a->t_csseed[0]/a->t_csseed[1]);
	fprintf(stderr, "supseed_time\t%.1f\t%.1f\t%.1fX\n",
		 a->t_supseed[0], a->t_supseed[1], a->t_supseed[0]/a->t_supseed[1]);
	fprintf(stderr, "input_reads\t%ld\n", a->reads_n);
	fprintf(stderr, "read_len\t%ld\n", a->read_len);
	fprintf(stderr, "Compression ration  \t%.2f\n", 1.0*a->reads_n*a->read_len/a->cs_bases);
	fprintf(stderr, "Mismatches number   \t%.2f\n", 1.0*a->mis_n[0]/a->reads_n);
	fprintf(stderr, "Reused seeds        \t%.2f\n", 1.0*a->reused_n[0]/a->reads_n);
	fprintf(stderr, "Supplementary seeds \t%.2f\n", 1.0*a->sup_n[0]/a->reads_n);
}

#endif //ZIP_SEEDING_5_9_TIME_PROF_H
