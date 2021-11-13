//
// Created by ixiaohu on 2021/10/14.
//

#include <stdio.h>

#include "cstl/kvec.h"
#include "cstl/kthread.h"
#include "bwalib/bwa.h"
#include "bwalib/kseq.h"
#include "bwalib/utils.h"
#include "samop.h"

KSEQ_DECLARE(gzFile)

typedef struct {
	int input_reads_bound;
	FILE *zip_seeds, *mem_seeds;
	gzFile zip_sam, mem_sam;
	int64_t n_processed;
} ktp_aux_t;

typedef struct {
	ktp_aux_t *aux;
	int n_seqs;
	uint8_t **zip_seeds, **mem_seeds;
	sam_line_t *zip_sam, *mem_sam;
} ktp_data_t;

typedef struct {
	long mem_seed_n, hit_seed_n;
	long mem_align_q3, hit_p_q3, hit_nm_q3;
	long mem_align_q10, hit_p_q10, hit_nm_q10;
	long mem_align_q20, hit_p_q20, hit_nm_q20;
	long mem_align_q30, hit_p_q30, hit_nm_q30;
	long mem_align_q40, hit_p_q40, hit_nm_q40;
} prof_t;
prof_t prof[256];

static uint8_t** load_seeds(int n_chunk_bound, FILE *f, int *n, long *size) {
	int i;
	*n = n_chunk_bound;
	*size = 0;
	uint8_t **seeds = malloc(n_chunk_bound * sizeof(uint8_t*));
	for (i = 0; i < n_chunk_bound; i++) {
		long bytes;
		fread(&bytes, sizeof(long), 1, f);
		if (feof(f)) {
			*n = i;
			break;
		}
		*size += bytes;
		seeds[i] = malloc(bytes * sizeof(uint8_t) + sizeof(long));
		*(long*)seeds[i] = bytes;
		fread(seeds[i] + sizeof(long), sizeof(uint8_t), bytes, f);
	}
	return seeds;
}

typedef struct {
	int qb, qe;
	long rb, re;
	long sa_size;
} seed_t;
typedef kvec_t(seed_t) seed_v;

static void seed_consistency(const uint8_t *mem_buf, const uint8_t *zip_buf, int t_id) {
	long bytes, p;
	int i, j, n;

	seed_v *mem_seeds = calloc(1, sizeof(seed_v));
	mem_seeds->n = 0;
	p = 0;
	bytes = *(long*)(mem_buf + p); p += sizeof(long);
	n = *(int*)(mem_buf + p); p += sizeof(int);
	for (i = 0; i < n; i++) {
		int qb = *(int*)(mem_buf + p); p += sizeof(int);
		int l = *(int*)(mem_buf + p); p += sizeof(int);
		int qe = qb + l;
		long sa_size = *(long*)(mem_buf + p); p += sizeof(long);
		int occ = sa_size < 500 ? sa_size : 500;
		for (j = 0; j < occ; j++) { // reference occurrence of mem seeds
			long rb = *(long*)(mem_buf + p); p += sizeof(long);
			long re = rb + l;
			seed_t s;
			s.qb = qb; s.qe = qe;
			s.rb = rb; s.re = re;
			s.sa_size = sa_size;
			kv_push(seed_t, *mem_seeds, s);
		}
	}
	assert(p == bytes + sizeof(long));

	seed_v *zip_seeds = calloc(1, sizeof(seed_v));
	zip_seeds->n = 0;
	p = 0;
	bytes = *(long*)(zip_buf + p); p += sizeof(long);
	n = *(int*)(zip_buf + p); p += sizeof(int);
	for (i = 0; i < n; i++) {
		int qb = *(int*)(zip_buf + p); p += sizeof(int);
		int l = *(int*)(zip_buf + p); p += sizeof(int);
		int qe = qb + l;
		long sa_size = *(long*)(zip_buf + p); p += sizeof(long);
		int occ = sa_size < 500 ? sa_size : 500;
		for (j = 0; j < occ; j++) { // reference occurrence of mem seeds
			long rb = *(long*)(zip_buf + p); p += sizeof(long);
			long re = rb + l;
			seed_t s;
			s.qb = qb; s.qe = qe;
			s.rb = rb; s.re = re;
			s.sa_size = sa_size;
			kv_push(seed_t, *zip_seeds, s);
		}
	}
	assert(p == bytes + sizeof(long));

	int non_rep_seeds = 0, cnt = 0;
	for (i = 0; i < mem_seeds->n; i++) {
		const seed_t *m = &mem_seeds->a[i];
		if (m->sa_size > 100) continue; // Ignore the repetitive seeds
		non_rep_seeds++;
		for (j = 0; j < zip_seeds->n; j++) {
			const seed_t *t = &zip_seeds->a[j];
			if (t->sa_size > 150) continue;
			if (t->qe <= m->qb || t->qb >= m->qe) continue;
			int que_ol = (t->qe < m->qe ?t->qe :m->qe) - (t->qb > m->qb ?t->qb :m->qb);
			int ref_ol = (t->re < m->re ?t->re :m->re) - (t->rb > m->rb ?t->rb :m->rb);
			int min_l  = (t->qe - t->qb < m->qe - m->qb) ?(t->qe - t->qb) :(m->qe - m->qb);
			if (que_ol >= (int)(min_l * 0.80 + 0.499) && ref_ol >= (int)(min_l * 0.80 + 0.499)) {
				cnt++;
				break;
			}
		}
	}
	prof[t_id].mem_seed_n += non_rep_seeds;
	prof[t_id].hit_seed_n += cnt;
	free(zip_seeds->a); free(zip_seeds);
	free(mem_seeds->a); free(mem_seeds);
}

static void worker1(void *_data, long seq_id, int t_id) {
	ktp_data_t *data = (ktp_data_t*)_data;
	seed_consistency(data->mem_seeds[seq_id], data->zip_seeds[seq_id], t_id);
	free(data->mem_seeds[seq_id]);
	free(data->zip_seeds[seq_id]);
}

static void *tp_check_seeds(void *_aux, int step, void *_data) {
	ktp_aux_t *aux = (ktp_aux_t*)_aux;
	if (step == 0) {
		double rtime_s = realtime();
		int n1, n2;
		long zip_size = 0, mem_size = 0;
		ktp_data_t *ret = calloc(1, sizeof(ktp_data_t));;
		ret->zip_seeds = load_seeds(aux->input_reads_bound, aux->zip_seeds, &n1, &zip_size);
		ret->mem_seeds = load_seeds(aux->input_reads_bound, aux->mem_seeds, &n2, &mem_size);
		assert(n1 == n2);
		ret->n_seqs = n1;
		if (ret->n_seqs == 0) {
			free(ret->zip_seeds);
			free(ret->mem_seeds);
			free(ret);
			return 0;
		}
		double rtime_e = realtime();
		fprintf(stderr, "[%s_step1] Input %ldM zip seeds %ldM mem seeds in %.1f real sec\n",
		  __func__, zip_size/1000/1000, mem_size/1000/1000, rtime_e-rtime_s);
		return ret;
	} else if (step == 1) {
		double rtime_s = realtime();
		ktp_data_t *data = (ktp_data_t*)_data;
		kt_for(16, worker1, data, data->n_seqs);
		free(data->zip_seeds); free(data->mem_seeds);
		aux->n_processed += data->n_seqs;
		double rtime_e = realtime();
		fprintf(stderr, "[%s_step2] Check seeds for %d reads in %.1f real sec, have processed %ld reads in total\n",
		        __func__, data->n_seqs, rtime_e-rtime_s, aux->n_processed);
		free(data);
		return 0;
	}
	return 0;
}

void check_seeds(const char *zip_fn, const char *mem_fn) {
	ktp_aux_t aux; memset(&aux, 0, sizeof(aux));
	aux.zip_seeds = fopen(zip_fn, "rb"); assert(aux.zip_seeds != NULL);
	aux.mem_seeds = fopen(mem_fn, "rb"); assert(aux.mem_seeds != NULL);
	aux.input_reads_bound = 1600 * 1000; // 160M bases for reads length of 100
	memset(prof, 0, sizeof(prof));
	kt_pipeline(2, tp_check_seeds, &aux, 2);
	int i;
	for (i = 1; i < 256; i++) prof[0].mem_seed_n += prof[i].mem_seed_n;
	for (i = 1; i < 256; i++) prof[0].hit_seed_n += prof[i].hit_seed_n;
	const prof_t *p = &prof[0];
	fprintf(stderr, "Seeds consistency: %.4f %% = %ld / %ld\n", 100.0 * p->hit_seed_n / p->mem_seed_n, p->hit_seed_n, p->mem_seed_n);
	fclose(aux.zip_seeds); fclose(aux.mem_seeds);
}

/************************
 * Check for alignment  *
 ************************/

static sam_line_t* load_sam(int n_chunk_bound, gzFile f, int *n) {
	int i;
	*n = n_chunk_bound;
	sam_line_t *sam = malloc(n_chunk_bound * sizeof(sam_line_t));
	for (i = 0; i < n_chunk_bound; i++) {
		sam_line_t zip_sam = fetch_samline1(f);
		if (gzeof(f)) {
			*n = i;
			break;
		}
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
	if (*n == 0) {
		free(sam);
		return 0;
	}
	return sam;
}

static void* tp_check_sam(void *_aux, int step, void *_data) {
	ktp_aux_t *aux = (ktp_aux_t*)_aux;
	if (step == 0) {
		double rtime_s = realtime();
		ktp_data_t *data = calloc(1, sizeof(ktp_data_t));
		int n1, n2;
		data->zip_sam = load_sam(aux->input_reads_bound, aux->zip_sam, &n1);
		data->mem_sam = load_sam(aux->input_reads_bound, aux->mem_sam, &n2);
		assert(n1 == n2);
		data->n_seqs = n1;
		if (data->n_seqs == 0) {
			free(data);
			return 0;
		}
		double rtime_e = realtime();
		fprintf(stderr, "[%s_step1] Input %d bwamem and zipmem primary alignments in %.1f real sec\n",
		        __func__, data->n_seqs, rtime_e-rtime_s);
		return data;
	} else if (step == 1) {
		double rtime_s = realtime();
		int i;
		ktp_data_t *data = (ktp_data_t*)_data;
		for (i = 0; i < data->n_seqs; i++) {
			const sam_line_t *m = &data->mem_sam[i];
			const sam_line_t *z = &data->zip_sam[i];
			assert(!strcmp(m->qname, z->qname));
			prof_t *p = &prof[0];
			int pos_true = (!strcmp(z->rname, m->rname)) && (abs(z->pos - m->pos) <= 5);
			int nm_true = z->nm == m->nm;
			if (m->mapq > 3) {
				p->mem_align_q3++;
				p->hit_p_q3 += pos_true;
				p->hit_nm_q3 += nm_true;
			}
			if (m->mapq > 10) {
				p->mem_align_q10++;
				p->hit_p_q10 += pos_true;
				p->hit_nm_q10 += nm_true;
			}
			if (m->mapq > 20) {
				p->mem_align_q20++;
				p->hit_p_q20 += pos_true;
				p->hit_nm_q20 += nm_true;
			}
			if (m->mapq > 30) {
				p->mem_align_q30++;
				p->hit_p_q30 += pos_true;
				p->hit_nm_q30 += nm_true;
			}
			if (m->mapq > 40) {
				p->mem_align_q40++;
				p->hit_p_q40 += pos_true;
				p->hit_nm_q40 += nm_true;
			}
			free(z->data); free(m->data);
		}
		free(data->zip_sam); free(data->mem_sam);
		aux->n_processed += data->n_seqs;
		double rtime_e = realtime();
		fprintf(stderr, "[%s_step2] Check alignments in %.1f real sec, have processed %ld reads in total\n",
		        __func__, rtime_e-rtime_s, aux->n_processed);
		return 0;
	}
	return 0;
}

void check_sam(const char *zip_fn, const char *mem_fn) {
	ktp_aux_t aux; memset(&aux, 0, sizeof(aux));
	aux.zip_sam = gzopen(zip_fn, "r"); assert(aux.zip_sam != NULL);
	aux.mem_sam = gzopen(mem_fn, "r"); assert(aux.mem_sam != NULL);
	aux.input_reads_bound = 1600 * 1000;
	memset(prof, 0, sizeof(prof));
	kt_pipeline(2, tp_check_sam, &aux, 2);
	gzclose(aux.zip_sam); gzclose(aux.mem_sam);

	const prof_t *p = &prof[0];
	fprintf(stderr, "Primary for Q3:  %ld\n", p->mem_align_q3);
	fprintf(stderr, "    POS:         %ld\t%.4f [%%]\n", p->hit_p_q3, 100.0 * p->hit_p_q3 / p->mem_align_q3);
	fprintf(stderr, "    NM:          %ld\t%.4f [%%]\n", p->hit_nm_q3, 100.0 * p->hit_nm_q3 / p->mem_align_q3);
	fprintf(stderr, "Primary for Q10: %ld\n", p->mem_align_q10);
	fprintf(stderr, "    POS:         %ld\t%.4f [%%]\n", p->hit_p_q10, 100.0 * p->hit_p_q10 / p->mem_align_q10);
	fprintf(stderr, "    NM:          %ld\t%.4f [%%]\n", p->hit_nm_q10, 100.0 * p->hit_nm_q10 / p->mem_align_q10);
	fprintf(stderr, "Primary for Q20: %ld\n", p->mem_align_q20);
	fprintf(stderr, "    POS:         %ld\t%.4f [%%]\n", p->hit_p_q20, 100.0 * p->hit_p_q20 / p->mem_align_q20);
	fprintf(stderr, "    NM:          %ld\t%.4f [%%]\n", p->hit_nm_q20, 100.0 * p->hit_nm_q20 / p->mem_align_q20);
	fprintf(stderr, "Primary for Q30: %ld\n", p->mem_align_q30);
	fprintf(stderr, "    POS:         %ld\t%.4f [%%]\n", p->hit_p_q30, 100.0 * p->hit_p_q30 / p->mem_align_q30);
	fprintf(stderr, "    NM:          %ld\t%.4f [%%]\n", p->hit_nm_q30, 100.0 * p->hit_nm_q30 / p->mem_align_q30);
	fprintf(stderr, "Primary for Q40: %ld\n", p->mem_align_q40);
	fprintf(stderr, "    POS:         %ld\t%.4f [%%]\n", p->hit_p_q40, 100.0 * p->hit_p_q40 / p->mem_align_q40);
	fprintf(stderr, "    NM:          %ld\t%.4f [%%]\n", p->hit_nm_q40, 100.0 * p->hit_nm_q40 / p->mem_align_q40);
}

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Check consistency of seeds or SAM alignments between zipmem and bwamem\n");
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "    check seed <zipmem.seeds> <bwamem.seeds>\n");
	fprintf(stderr, "    check sam <zipmem.sam> <bwamem.sam>\n");
	fprintf(stderr, "\n");
	return 0;
}

int main(int argc, char *argv[]) {
	if (argc == 1) return usage();
	if (!strcmp(argv[1], "seed") && argc == 4) check_seeds(argv[2], argv[3]);
	else if (!strcmp(argv[1], "sam") && argc == 4) check_sam(argv[2], argv[3]);
	else return usage();
	return 0;
}
