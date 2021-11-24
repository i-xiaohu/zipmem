#include <stdio.h>
#include "cstl/kthread.h"
#include "bwalib/bwa.h"
#include "bwalib/kseq.h"
#include "bwalib/kopen.h"
#include "bwalib/utils.h"
#include "time_prof.h"
#include "zipmem.h"
#include "read_input.h"

#define ZIPMEM_VERSION "5.2-m23"

KSEQ_DECLARE(gzFile)

zsmem_prof_t zmp;

typedef struct {
	kseq_t *ks;
	long n_has_input;
	gzFile freads; // Only reading bases
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
		if (aux->ks) ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, NULL);
		else ret->seqs = load_reads(aux->n_has_input, aux->actual_chunk_size, &ret->n_seqs, aux->freads);
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		aux->n_has_input += ret->n_seqs;
		if (!aux->copy_comment)
			for (i = 0; i < ret->n_seqs; ++i) {
				free(ret->seqs[i].comment);
				ret->seqs[i].comment = 0;
			}
		for (i = 0; i < ret->n_seqs; ++i) size += ret->seqs[i].l_seq;
		return ret;
	} else if (step == 1) {
		double ctime_s = cputime(), rtime_s = realtime();
		const mem_opt_t *opt = aux->opt;
		const bwaidx_t *idx = aux->idx;
		zipmem_seeding(opt, idx->bwt, data->n_seqs, data->seqs);
		aux->n_processed += data->n_seqs;
		double ctime_e = cputime(), rtime_e = realtime();
		zmp.t_seeding[0] += ctime_e - ctime_s; zmp.t_seeding[1] += rtime_e - rtime_s;
		return data;
	} else if (step == 2) {
		for (i = 0; i < data->n_seqs; ++i) {
			// sam here is used to stores the seeds and the preceding number indicating memory cost.
			long bytes = *(long*)data->seqs[i].sam + sizeof(long);
			fwrite(data->seqs[i].sam, sizeof(uint8_t), bytes, stdout);
			free(data->seqs[i].sam);
			free(data->seqs[i].name); free(data->seqs[i].comment);
			free(data->seqs[i].seq); free(data->seqs[i].qual);
		}
		free(data->seqs); free(data);
		fprintf(stderr, "\033[K");
		fprintf(stderr, "\033[33mSeeding for %ld reads costs %.1f CPU seconds, %.1f real seconds.\033[0m\n",
	        aux->n_processed, zmp.t_seeding[0], zmp.t_seeding[1]);
		fprintf(stderr, "\033[1A");
		return 0;
	}
	return 0;
}

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Zipmem (version %s) is the compressive version of `bwamem seeding` for reordered reads\n", ZIPMEM_VERSION);
	fprintf(stderr, "Usage: zipmem [options] <FM-index> data.reads > zip.seeds\n");
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
	fprintf(stderr, "Use `bwamem extend` to extend the zipmem seeds into final alignments\n");
	fprintf(stderr, "\n");
	return 0;
}

int main(int argc, char *argv[]) {
	if (argc == 1) return usage();
	double rtime = realtime();
	int i, c;
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
		else if (c == 't') opt->n_threads = (int)strtol(optarg, NULL, 10), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'K') fixed_chunk_size = (int)strtol(optarg, NULL, 10);
		else if (c == 'v') bwa_verbose = (int)strtol(optarg, NULL, 10);
		else if (c == '1') no_mt_io = 1;
		else return 1;
	}

	if (optind + 1 >= argc || optind + 2 < argc) {
		fprintf(stderr, "[%s] Invalid usage\n", __func__ );
		free(opt);
		return 1;
	}

	aux.idx = bwa_idx_load_from_shm(argv[optind]);
	if (aux.idx == 0) {
		if ((aux.idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1;
	} else if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] load the bwa index from shared memory\n", __func__);

	int not_fastaq = is_reads_file(argv[optind + 1]);
	if (not_fastaq) {
		aux.freads = gzopen(argv[optind + 1], "r");
	} else {
		ko = kopen(argv[optind + 1], &fd);
		if (ko == 0) {
			if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
			return 1;
		}
		fp = gzdopen(fd, "r");
		aux.ks = kseq_init(fp);
	}

	memset(&zmp, 0, sizeof(zmp));
	aux.actual_chunk_size = fixed_chunk_size > 0? fixed_chunk_size : opt->chunk_size * opt->n_threads;
	kt_pipeline(no_mt_io? 1 : 2, tp_seeding, &aux, 3);
	fprintf(stderr, "\n");
//	zsmem_prof_output(&zmp);

	free(opt);
	bwa_idx_destroy(aux.idx);
	if (not_fastaq) {
		gzclose(aux.freads);
	} else {
		kseq_destroy(aux.ks);
		err_gzclose(fp); kclose(ko);
	}

	err_fflush(stdout);
	err_fclose(stdout);
	fprintf(stderr, "Zipmem Version: %s\n", ZIPMEM_VERSION);
	fprintf(stderr, "CMD:");
	for (i = 0; i < argc; ++i) fprintf(stderr, " %s", argv[i]); fprintf(stderr, "\n");
	fprintf(stderr, "Note: The time cost of writing seeds into disk is excluded. Because it does not happen in real alignment scenario, ");
	fprintf(stderr, "where seeds are buffered in memory.\n");
	return 0;
}