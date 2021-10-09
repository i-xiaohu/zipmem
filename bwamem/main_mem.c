#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "../bwalib/bwa.h"
#include "../bwalib/utils.h"
#include "../bwalib/kopen.h"
#include "../bwalib/kseq.h"
#include "../cstl/kthread.h"
#include "bwamem.h"

#ifndef BWA_PACKAGE_VERSION
#define BWA_PACKAGE_VERSION "0.7.17-r1188"
#endif

KSEQ_DECLARE(gzFile)

typedef struct {
	kseq_t *ks, *ks2;
	FILE *fseeds;
	mem_opt_t *opt;
	mem_pestat_t *pes0;
	int64_t n_processed;
	int copy_comment, actual_chunk_size;
	bwaidx_t *idx;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
	uint8_t **seeds;
} ktp_data_t;

// Three steps of pipeline for seeding: input reads, seeding, output seeds.
static void *tp_seeding(void *shared, int step, void *_data) {
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if (step == 0) {
		ktp_data_t *ret;
		int64_t size = 0;
		ret = calloc(1, sizeof(ktp_data_t));
		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2);
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		if (!aux->copy_comment)
			for (i = 0; i < ret->n_seqs; ++i) {
				free(ret->seqs[i].comment);
				ret->seqs[i].comment = 0;
			}
		for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[%s_step1] Input %d reads (%ld bp)\n", __func__, ret->n_seqs, (long)size);
		return ret;
	} else if (step == 1) {
		const mem_opt_t *opt = aux->opt;
		const bwaidx_t *idx = aux->idx;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[%s_step2] Process reads\n", __func__ );
		mem_seeding(opt, idx->bwt, data->n_seqs, data->seqs);
		aux->n_processed += data->n_seqs;
		return data;
	} else if (step == 2) {
		for (i = 0; i < data->n_seqs; ++i) {
			// sam here is used to stores the seeds and the preceding number indicating memory cost.
			long bytes = *(long*)data->seqs[i].sam + sizeof(long);
			fwrite(data->seqs[i].sam, sizeof(uint8_t), bytes, stdout);
			free(data->seqs[i].name); free(data->seqs[i].comment);
			free(data->seqs[i].seq); free(data->seqs[i].qual); free(data->seqs[i].sam);
		}
		free(data->seqs); free(data);
		if (bwa_verbose >= 3)
			fprintf(stderr, "[%s_step3] Output seeds. %ld reads have been processed\n", __func__ , aux->n_processed);
		return 0;
	}
	return 0;
}

static int seeding_usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: bwamem seeding [options] <FM-index> <in1.fq> [in2.fq] > mem.seeds\n");
	fprintf(stderr, "Seeding options:\n");
	fprintf(stderr, "    -k INT    minimum seed length [19]\n");
	fprintf(stderr, "    -r FLOAT  look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]\n");
	fprintf(stderr, "    -y INT    seed occurrence for the 3rd round seeding [20]\n");
	fprintf(stderr, "    -c INT    skip seeds with more than INT occurrences [500]\n");
	fprintf(stderr, "    -s INT    MEM SA size upper bound to trigger re-seeding [10]\n");
	fprintf(stderr, "Common options(also suitable for `extend`):\n");
	fprintf(stderr, "    -t INT    number of threads [1]\n");
	fprintf(stderr, "    -K INT    process INT input bases in each batch regardless of nThreads (for reproducibility) [10M]\n");
	fprintf(stderr, "    -v INT    verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
	fprintf(stderr, "    -1        Disable multiple I/O threads.\n");
	fprintf(stderr, "\n");
	return 0;
}

int seeding_main(int argc, char *argv[]) {
	if (argc == 1) return seeding_usage();
	int c;
	int no_mt_io = 0; // Multiple IO thread label
	int fixed_chunk_size = -1;
	int fd, fd2;
	void *ko = 0, *ko2 = 0;
	gzFile fp, fp2 = 0;
	mem_opt_t *opt;
	ktp_aux_t aux;

	memset(&aux, 0, sizeof(ktp_aux_t));
	aux.opt = opt = mem_opt_init();
	while ((c = getopt(argc, argv, "k:r:y:c:s:t:K:v:1")) >= 0) {
		// Seeding options
		if (c == 'k') opt->min_seed_len = (int)strtol(optarg, NULL, 10);
		else if (c == 'r') opt->split_factor = strtof(optarg, NULL);
		else if (c == 'y') opt->max_mem_intv = strtol(optarg, NULL, 10);
		else if (c == 'c') opt->max_occ = (int)strtol(optarg, NULL, 10);
		else if (c == 's') opt->split_width = (int)strtol(optarg, NULL, 10);
		// Common options
		if (c == 't') opt->n_threads = (int)strtol(optarg, NULL, 10), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'K') fixed_chunk_size = (int)strtol(optarg, NULL, 10);
		else if (c == 'v') bwa_verbose = (int)strtol(optarg, NULL, 10);
		else if (c == '1') no_mt_io = 1;
		else return 1;
	}

	if (optind + 1 >= argc || optind + 3 < argc) {
		free(opt);
		return 1;
	}

	aux.idx = bwa_idx_load_from_shm(argv[optind]);
	if (aux.idx == 0) {
		if ((aux.idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);

	ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
		return 1;
	}
	fp = gzdopen(fd, "r");
	aux.ks = kseq_init(fp);
	if (optind + 2 < argc) {
		ko2 = kopen(argv[optind + 2], &fd2);
		if (ko2 == 0) {
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 2]);
			return 1;
		}
		fp2 = gzdopen(fd2, "r");
		aux.ks2 = kseq_init(fp2);
		opt->flag |= MEM_F_PE;
	}
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
	kt_pipeline(no_mt_io? 1 : 2, tp_seeding, &aux, 3);
	free(opt);
	bwa_idx_destroy(aux.idx);
	kseq_destroy(aux.ks);
	err_gzclose(fp); kclose(ko);
	if (aux.ks2) {
		kseq_destroy(aux.ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	return 0;
}

void output_seeds(int n, uint8_t **all_seeds, bseq1_t *reads) {
	int r, i, j;
	for (r = 0; r < n; r++) {
		long pointer = 0;
		uint8_t *seeds = all_seeds[r];
		int mem_n = *(int*)(seeds + pointer); pointer += sizeof(int);
		kstring_t tmp = {0, 0, 0};
		ksprintf(&tmp, "%d\n", mem_n);
		for (i = 0; i < mem_n; ++i) {
			int qb = *(int*)(seeds + pointer); pointer += sizeof(int);
			int slen = *(int*)(seeds + pointer); pointer += sizeof(int);
			long sa_size = *(long*)(seeds + pointer); pointer += sizeof(long);
			ksprintf(&tmp, "%d %d %ld\n", qb, slen, sa_size);
			int true_occ = sa_size < 500 ?sa_size :500;
			for (j = 0; j < true_occ; j++) {
				long pos = *(long*)(seeds + pointer); pointer += sizeof(long);
				ksprintf(&tmp, "%ld ", pos);
			}
			ksprintf(&tmp, "\n");
		}
		reads[r].sam = tmp.s;
	}
}

// Three steps of pipeline for extending: input seeds and reads, extending, output SAM.
static void *tp_extending(void *shared, int step, void *_data) {
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if (step == 0) {
		ktp_data_t *ret;
		int64_t size = 0;
		ret = calloc(1, sizeof(ktp_data_t));
		ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, aux->ks2);
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		if (!aux->copy_comment)
			for (i = 0; i < ret->n_seqs; ++i) {
				free(ret->seqs[i].comment);
				ret->seqs[i].comment = 0;
			}
		for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
		ret->seeds = malloc(sizeof(uint8_t*) * ret->n_seqs);
		long total_bytes = 0;
		for (i = 0; i < ret->n_seqs; i++) {
			long seed_bytes = 0;
			fread(&seed_bytes, sizeof(long), 1, aux->fseeds);
			total_bytes += seed_bytes;
			ret->seeds[i] = malloc(seed_bytes);
			fread(ret->seeds[i], sizeof(uint8_t), seed_bytes, aux->fseeds);
		}
		if (bwa_verbose >= 3) {
			fprintf(stderr, "[%s_step1] Input %d reads (%ld bp) and %ld MB seeds\n",
		            __func__, ret->n_seqs, (long)size, total_bytes / 1024 / 1024);
		}
		return ret;
	} else if (step == 1) {
		const mem_opt_t *opt = aux->opt;
		const bwaidx_t *idx = aux->idx;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[%s_step2] Process reads\n", __func__ );
//		output_seeds(data->n_seqs, data->seeds, data->seqs);
		mem_extend(opt, idx->bwt, idx->bns, idx->pac, aux->n_processed, data->n_seqs, data->seeds, data->seqs, aux->pes0);
		aux->n_processed += data->n_seqs;
		return data;
	} else if (step == 2) {
		for (i = 0; i < data->n_seqs; ++i) {
			if (data->seqs[i].sam) err_fputs(data->seqs[i].sam, stdout);
			free(data->seqs[i].name); free(data->seqs[i].comment);
			free(data->seqs[i].seq); free(data->seqs[i].qual);
			free(data->seqs[i].sam);
			// Seeds are deallocated by the core function in step2.
		}
		free(data->seqs); free(data);
		if (bwa_verbose >= 3)
			fprintf(stderr, "[%s_step3] Output SAM. %ld reads have been processed\n", __func__ , aux->n_processed);
		return 0;
	}
	return 0;
}

static void update_a(mem_opt_t *opt, const mem_opt_t *opt0) {
	if (opt0->a) { // matching score is changed
		if (!opt0->b) opt->b *= opt->a;
		if (!opt0->T) opt->T *= opt->a;
		if (!opt0->o_del) opt->o_del *= opt->a;
		if (!opt0->e_del) opt->e_del *= opt->a;
		if (!opt0->o_ins) opt->o_ins *= opt->a;
		if (!opt0->e_ins) opt->e_ins *= opt->a;
		if (!opt0->zdrop) opt->zdrop *= opt->a;
		if (!opt0->pen_clip5) opt->pen_clip5 *= opt->a;
		if (!opt0->pen_clip3) opt->pen_clip3 *= opt->a;
		if (!opt0->pen_unpaired) opt->pen_unpaired *= opt->a;
	}
}

static int extend_usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: bwamem extend [options] <FM-index> <mem.seeds> <in1.fq> [in2.fq] > aln.sam\n");
	fprintf(stderr, "Seeds chaining options\n");
	fprintf(stderr, "    -c INT    skip seeds with more than INT occurrences [500]\n");
	fprintf(stderr, "    -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.5]\n");
	fprintf(stderr, "    -W INT        discard a chain if seeded bases shorter than INT [0]\n");
	fprintf(stderr, "    -G INT        maximum gap between seeds to chain [10000]\n");
	fprintf(stderr, "    -N INT        maximum chains to extend [INF]\n");
	fprintf(stderr, "    -X FLOAT      mask level for significantly overlapped chain [0.5]\n");
	fprintf(stderr, "Smith-Waterman options:\n");
	fprintf(stderr, "    -w INT        band width for banded alignment [100]\n");
	fprintf(stderr, "    -d INT        off-diagonal X-dropoff [100]\n");
	fprintf(stderr, "    -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]\n");
	fprintf(stderr, "    -B INT        penalty for a mismatch [4]\n");
	fprintf(stderr, "    -O INT[,INT]  gap open penalties for deletions and insertions [6,6]\n");
	fprintf(stderr, "    -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]\n");
	fprintf(stderr, "    -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]\n");
	fprintf(stderr, "    -U INT        penalty for an unpaired read pair [17]\n");
	fprintf(stderr, "Paired-end options:\n");
	fprintf(stderr, "    -I FLOAT[,FLOAT[,INT[,INT]]]\n");
	fprintf(stderr, "                  specify the mean, standard deviation, max and min of the insert size distribution.\n");
	fprintf(stderr, "                  FR orientation only inferred.\n");
	fprintf(stderr, "    -m INT        perform at most INT rounds of mate rescues for each read [50]\n");
	fprintf(stderr, "    -S            skip mate rescue\n");
	fprintf(stderr, "    -P            skip pairing; mate rescue performed unless -S also in use\n");
	fprintf(stderr, "Formatting SAM options:\n");
	fprintf(stderr, "    -T INT        minimum score to output [30]\n");
	fprintf(stderr, "    -a            output all alignments for SE or unpaired PE\n");
	fprintf(stderr, "    -C            append FASTA/FASTQ comment to SAM output\n");
	fprintf(stderr, "\n");
	return 0;
}

int extend_main(int argc, char *argv[]) {
	if (argc == 1) return extend_usage();

	mem_opt_t *opt, opt0;
	int i, c, no_mt_io = 0;
	int fixed_chunk_size = -1;
	char *p;
	int fd, fd2;
	void *ko = 0, *ko2 = 0;
	gzFile fp, fp2 = 0;
	mem_pestat_t pes[4];
	ktp_aux_t aux;

	memset(&aux, 0, sizeof(ktp_aux_t));
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	aux.opt = opt = mem_opt_init();
	memset(&opt0, 0, sizeof(mem_opt_t));
	while ((c = getopt(argc, argv, "t:K:v:1c:D:W:G:N:X:w:d:A:B:O:E:L:U:I:m:SPT:aC")) >= 0) {
		if (c == 't') opt->n_threads = (int)strtol(optarg, NULL, 10), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'K') fixed_chunk_size = (int)strtol(optarg, NULL, 10);
		else if (c == 'v') bwa_verbose = (int)strtol(optarg, NULL, 10);
		else if (c == '1') no_mt_io = 1;
		// Chaining
		else if (c == 'c') opt->max_occ = (int)strtol(optarg, NULL, 10);
		else if (c == 'D') opt->drop_ratio = strtof(optarg, NULL);
		else if (c == 'W') opt->min_chain_weight = (int)strtol(optarg, NULL, 10);
		else if (c == 'G') opt->max_chain_gap = (int)strtol(optarg, NULL, 10);
		else if (c == 'N') opt->max_chain_extend = (int)strtol(optarg, NULL, 10);
		else if (c == 'X') opt->mask_level = strtof(optarg, NULL);
		// Smith-Waterman
		else if (c == 'w') opt->w = (int)strtol(optarg, NULL, 10);
		else if (c == 'd') opt->zdrop = (int)strtol(optarg, NULL, 10);
		else if (c == 'A') opt->a = (int)strtol(optarg, NULL, 10), opt0.a = 1;
		else if (c == 'B') opt->b = (int)strtol(optarg, NULL, 10), opt0.b = 1;
		else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = (int)strtol(optarg, &p, 10);
			if (*p != 0) opt->o_ins = (int)strtol(p+1, &p, 10);
		} else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = (int)strtol(optarg, &p, 10);
			if (*p != 0) opt->e_ins = (int)strtol(p+1, &p, 10);
		} else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = (int)strtol(optarg, &p, 10);
			if (*p != 0) opt->pen_clip3 = (int)strtol(p+1, &p, 10);
		}
		else if (c == 'U') opt->pen_unpaired = (int)strtol(optarg, &p, 10), opt0.pen_unpaired = 1;
		// Paired-end
		else if (c == 'I') { // specify the insert size distribution
			aux.pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0) pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0) pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0) pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
				        __func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		else if (c == 'm') opt->max_matesw = (int)strtol(optarg, &p, 10);
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		// SAM format
		else if (c == 'T') opt->T = (int)strtol(optarg, &p, 10), opt0.T = 1;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'C') aux.copy_comment = 1;
		else return 1;
	}

	if (optind + 1 >= argc || optind + 4 < argc) {
		free(opt);
		return 1;
	}

	update_a(opt, &opt0);
	bwa_fill_scmat(opt->a, opt->b, opt->mat);

	aux.idx = bwa_idx_load_from_shm(argv[optind]);
	if (aux.idx == 0) {
		if ((aux.idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);

	aux.fseeds = fopen(argv[optind + 1], "rb");
	ko = kopen(argv[optind + 2], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 2]);
		return 1;
	}
	fp = gzdopen(fd, "r");
	aux.ks = kseq_init(fp);
	if (optind + 3 < argc) {
		ko2 = kopen(argv[optind + 3], &fd2);
		if (ko2 == 0) {
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 3]);
			return 1;
		}
		fp2 = gzdopen(fd2, "r");
		aux.ks2 = kseq_init(fp2);
		opt->flag |= MEM_F_PE;
	}

	bwa_print_sam_hdr(aux.idx->bns, NULL);
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
	kt_pipeline(no_mt_io? 1 : 2, tp_extending, &aux, 3);
	free(opt);
	bwa_idx_destroy(aux.idx);
	fclose(aux.fseeds);
	kseq_destroy(aux.ks);
	err_gzclose(fp); kclose(ko);
	if (aux.ks2) {
		kseq_destroy(aux.ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	return 0;
}

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: bwamem <seeding / extend>\n");
	fprintf(stderr, "  BWA-MEM (version %s) follows the seed-and-extend paradigm\n", BWA_PACKAGE_VERSION);
	fprintf(stderr, "  Command `seeding` is collecting candidate hits for query reads\n");
	fprintf(stderr, "  Command `extend` performs Smith-Waterman around seeds for final alignments\n");
	fprintf(stderr, "  In BWA source, they are an entity, we split them for testing zip-seeding\n");
	fprintf(stderr, "  It does not impact the original performance of BWA-MEM, users can check the time cost or SAM files\n");
	fprintf(stderr, "  Sample usage:\n");
	fprintf(stderr, "    bwamem seeding <FM-index> data1.fq > data1.seeds\n");
	fprintf(stderr, "    bwamem extend <FM-index> data1.seeds data1.fq > data.sam\n");
	fprintf(stderr, "\n");
	return 0;
}

int main(int argc, char *argv[]) {
	if (argc == 1) return usage();
	int i, ret;
	extern char *bwa_pg;
	kstring_t pg = {0,0,0};
	ksprintf(&pg, "@PG\tID:bwamem\tPN:bwamem\tVN:%s\tCL:%s", BWA_PACKAGE_VERSION, argv[0]);
	for (i = 1; i < argc; ++i) ksprintf(&pg, " %s", argv[i]);
	bwa_pg = pg.s;

	double rtime = realtime();
	if (!strcmp(argv[1], "seeding")) ret = seeding_main(argc-1, argv+1);
	else if (!strcmp(argv[1], "extend")) ret = extend_main(argc-1, argv+1);
	else {
		ret = 1;
		fprintf(stderr, "Only commands `seeding` and `extend` are allowed\n");
	}
	err_fflush(stdout);
	err_fclose(stdout);
	if (ret == 0) {
		fprintf(stderr, "BWA-MEM Version: %s\n", BWA_PACKAGE_VERSION);
		fprintf(stderr, "CMD:");
		for (i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n");
		fprintf(stderr, "Time cost for %s: %.2f real sec, %.2f CPU sec\n", argv[1], realtime()-rtime, cputime());
	}
	free(bwa_pg);
	return 0;
}