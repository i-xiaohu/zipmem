//
// Created by ixiaohu on 2021/11/13.
//

#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>

#include "bwalib/bwa.h"
#include "bwalib/kseq.h"
#include "bwalib/utils.h"
#include "bwalib/kopen.h"
#include "cstl/kthread.h"
#include "cstl/kvec.h"
#include "samop.h"
#include "read_input.h"

KSEQ_DECLARE(gzFile)

typedef struct {
	int actual_chunk_size;
	kseq_t *ks;
	long n_has_input;
	gzFile freads; // Only reading bases
	FILE *zip_seeds, *mem_seeds;
	gzFile zip_sam, mem_sam;
	int64_t n_processed;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	bseq1_t *seqs;
	uint8_t **zip_seeds, **mem_seeds;
	sam_line_t *zip_sam, *mem_sam;
} ktp_data_t;

static uint8_t** load_seeds(int n, FILE *f) {
	int i;
	uint8_t **seeds = malloc(n * sizeof(uint8_t*));
	for (i = 0; i < n; i++) {
		long bytes;
		fread(&bytes, sizeof(long), 1, f);
		assert(!feof(f));
		seeds[i] = malloc(bytes * sizeof(uint8_t) + sizeof(long));
		*(long*)seeds[i] = bytes;
		fread(seeds[i] + sizeof(long), sizeof(uint8_t), bytes, f);
	}
	return seeds;
}

static sam_line_t* load_sam(int n, gzFile f) {
	int i;
	sam_line_t *sam = malloc(n * sizeof(sam_line_t));
	for (i = 0; i < n; i++) {
		sam_line_t zip_sam = fetch_samline1(f);
		assert(!gzeof(f));
		if (zip_sam.data == NULL) {
			i--;
			continue;
		}
		if (sec_ali(zip_sam.flag) || sup_ali(zip_sam.flag)) {
			free(zip_sam.data);
			i--;
			continue;
		}
		sam[i] = zip_sam;
	}
	return sam;
}

typedef struct {
	int qb, qe;
} seed_t;

typedef kvec_t(seed_t) seed_v;

static seed_v parse_encoded_seeds(const uint8_t *enc_seeds) {
	seed_v seeds; kv_init(seeds);
	long p = 0;
	long bytes = *(long*)(enc_seeds + p); p += sizeof(long);
	int i, j, n = *(int*)(enc_seeds + p); p += sizeof(int);
	for (i = 0; i < n; i++) {
		int qb = *(int*)(enc_seeds + p); p += sizeof(int);
		int l = *(int*)(enc_seeds + p); p += sizeof(int);
		int qe = qb + l;
		seed_t s; s.qb = qb; s.qe = qe; kv_push(seed_t, seeds, s);
		long sa_size = *(long*)(enc_seeds + p); p += sizeof(long);
		int occ = sa_size < 500 ? sa_size : 500;
		p += occ * sizeof(long);
	}
	assert(p == bytes + sizeof(long));
	return seeds;
}

static void process(const uint8_t *zseed, const uint8_t *mseed, const sam_line_t *zsam, const sam_line_t *msam) {
	int i;
	if (zsam->pos != msam->pos && msam->mapq > 40) {
		seed_v mem_seeds = parse_encoded_seeds(mseed);
		seed_v zip_seeds = parse_encoded_seeds(zseed);
		fprintf(stdout, "Read: %s\n", zsam->qname); assert(!strcmp(msam->qname, zsam->qname));
		fprintf(stdout, "Seq:  %s\n", zsam->seq);
		fprintf(stdout, "RNAME(ZM): %s\t%s\n", zsam->rname, msam->rname);
		fprintf(stdout, "FLAG(ZM):  %d\t%d\n", zsam->flag, msam->flag);
		fprintf(stdout, "POS(ZM):   %d\t%d\n", zsam->pos, msam->pos);
		fprintf(stdout, "MAPQ(ZM):  %d\t%d\n", zsam->mapq, msam->mapq);
		fprintf(stdout, "AS(ZM):    %d\t%d\n", zsam->as, msam->as);
		fprintf(stdout, "NM(ZM):    %d\t%d\n", zsam->nm, msam->nm);
		fprintf(stdout, "CIGAR(ZM): %s\t%s\n", zsam->cigar, msam->cigar);
		fprintf(stdout, "Seeds(Z):  ");
		for (i = 0; i < zip_seeds.n; i++)
			fprintf(stdout, "[%d, %d)\t", zip_seeds.a[i].qb, zip_seeds.a[i].qe); fprintf(stdout, "\n");
		fprintf(stdout, "Seeds(M):  ");
		for (i = 0; i < mem_seeds.n; i++)
			fprintf(stdout, "[%d, %d)\t", mem_seeds.a[i].qb, mem_seeds.a[i].qe); fprintf(stdout, "\n");
		fprintf(stdout, "\n");
	}
}

// Three steps of pipeline for seeding: input reads, seeding, output seeds.
static void *tp_study(void *shared, int step, void *_data) {
	ktp_aux_t *aux = (ktp_aux_t*)shared;
	ktp_data_t *data = (ktp_data_t*)_data;
	int i;
	if (step == 0) {
		double rtime_s = realtime();
		ktp_data_t *ret;
		int64_t size = 0;
		ret = calloc(1, sizeof(ktp_data_t));
		if (aux->ks) ret->seqs = bseq_read(aux->actual_chunk_size, &ret->n_seqs, aux->ks, NULL);
		else ret->seqs = load_reads(aux->n_has_input, aux->actual_chunk_size, &ret->n_seqs, aux->freads);
		if (ret->seqs == 0) {
			free(ret);
			return 0;
		}
		ret->mem_seeds = load_seeds(ret->n_seqs, aux->mem_seeds);
		ret->zip_seeds = load_seeds(ret->n_seqs, aux->zip_seeds);
		ret->mem_sam = load_sam(ret->n_seqs, aux->mem_sam);
		ret->zip_sam = load_sam(ret->n_seqs, aux->zip_sam);
		aux->n_has_input += ret->n_seqs;
		for (i = 0; i < ret->n_seqs; ++i) {
			free(ret->seqs[i].comment);
			ret->seqs[i].comment = 0;
			size += ret->seqs[i].l_seq;
		}
		double rtime_e = realtime();
		fprintf(stderr, "[%s_step1] Input done in %.1f real sec\n", __func__, rtime_e-rtime_s);
		return ret;
	} else if (step == 1) {
		double rtime_s = realtime();
		for (i = 0; i < data->n_seqs; i++) {
			process(data->zip_seeds[i], data->mem_seeds[i], &data->zip_sam[i], &data->mem_sam[i]);
		}
		aux->n_processed += data->n_seqs;
		double rtime_e = realtime();
		fprintf(stderr, "[%s_step2] Study done in %.1f real sec\n", __func__,  rtime_e-rtime_s);
		return data;
	} else if (step == 2) {
		for (i = 0; i < data->n_seqs; ++i) {
			free(data->seqs[i].name); free(data->seqs[i].comment);
			free(data->seqs[i].seq); free(data->seqs[i].qual);
			free(data->zip_seeds[i]); free(data->mem_seeds[i]);
			free(data->zip_sam[i].data);  free(data->mem_sam[i].data);
		}
		free(data->seqs);
		free(data->zip_seeds); free(data->mem_seeds);
		free(data->zip_sam); free(data->mem_sam);
		free(data);
		fprintf(stderr, "[%s_step3] %ld reads have been processed\n", __func__ , aux->n_processed);
		return 0;
	}
	return 0;
}

int main(int argc, char *argv[]) {
	if (argc != 6) {
		fprintf(stderr, "Usage: study_diff <reads> <mem_seeds> <zip_seeds> <mem_sam> <zip_sam>\n");
		return 0;
	}
	ktp_aux_t aux; memset(&aux, 0, sizeof(aux));

	int fd, not_fastaq = is_reads_file(argv[1]);
	gzFile fp; void *ko;
	if (not_fastaq) {
		aux.freads = gzopen(argv[1], "r");
	} else {
		ko = kopen(argv[1], &fd);
		fp = gzdopen(fd, "r");
		aux.ks = kseq_init(fp);
	}
	aux.mem_seeds = fopen(argv[2], "r");
	aux.zip_seeds = fopen(argv[3], "r");
	aux.mem_sam = gzopen(argv[4], "r");
	aux.zip_sam = gzopen(argv[5], "r");
	aux.actual_chunk_size = 120 * 1000 * 1000;

	kt_pipeline(2, tp_study, &aux, 3);

	if (not_fastaq) {
		gzclose(aux.freads);
	} else {
		kseq_destroy(aux.ks);
		err_gzclose(fp); kclose(ko);
	}
	return 0;
}