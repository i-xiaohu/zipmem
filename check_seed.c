//
// Created by ixiaohu on 2021/10/14.
//

#include <stdio.h>

#include "cstl/kvec.h"
#include "cstl/kthread.h"
#include "bwalib/bwa.h"
#include "bwalib/kseq.h"
#include "bwalib/utils.h"

KSEQ_DECLARE(gzFile)

typedef struct {
	char *qname;
	int flag;
	char *rname;
	int pos;
	int mapq;
	char *cigar;
	char *rnext;
	int pnext;  // position of next read
	int tlen;   // inferred insert size
	char *seq;
	char *qual;
	// Options
	int nm;
	int as;
	char *data; // All fields are allocated in one memory chunk.
} sam_line_t;

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
	long mem_align_n, hit_p_n, hit_nm_n, hit_mq_n; // Only primary alignments with MAPQ>3 of BWA-MEM are considered
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
		if (m->sa_size > 150) continue; // Ignore the repetitive seeds
		non_rep_seeds++;
		for (j = 0; j < zip_seeds->n; j++) {
			const seed_t *t = &zip_seeds->a[j];
			if (t->sa_size > 150) continue;
			if (t->qe <= m->qb || t->qb >= m->qe) continue;
			int que_ol = (t->qe < m->qe ?t->qe :m->qe) - (t->qb > m->qb ?t->qb :m->qb);
			int ref_ol = (t->re < m->re ?t->re :m->re) - (t->rb > m->rb ?t->rb :m->rb);
			int min_l  = (t->qe - t->qb < m->qe - m->qb) ?(t->qe - t->qb) :(m->qe - m->qb);
			if (que_ol >= (int)(min_l * 0.75 + 0.499) && ref_ol >= (int)(min_l * 0.75 + 0.499)) {
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
		fprintf(stderr, "[%s_step1] Input %ld bytes zip seeds %ld bytes mem seeds in %.1f real sec\n",
		  __func__, zip_size, mem_size, rtime_e-rtime_s);
		return ret;
	} else if (step == 1) {
		double rtime_s = realtime();
		ktp_data_t *data = (ktp_data_t*)_data;
		int i, j;
		long mem_substring = 0, zip_substring = 0;
		long mem_occ = 0, zip_occ = 0;
		for (i = 0; i < data->n_seqs; i++) {
			uint8_t *m = data->mem_seeds[i];
			long p = sizeof(long);
			int mem_n = *(int*)(m + p); p += sizeof(int);
			mem_substring += mem_n;
			for (j = 0; j < mem_n; j++) {
				p += sizeof(int);
				p += sizeof(int);
				long sa_size = *(long*)(m + p); p += sizeof(long);
				int occ = sa_size < 500 ? sa_size : 500;
				p += occ * sizeof(long);
				mem_occ += occ;
			}
			free(m);

			uint8_t *z = data->zip_seeds[i];
			p = sizeof(long);
			int zip_n = *(int*)(z + p); p += sizeof(int);
			zip_substring += zip_n;
			for (j = 0; j < zip_n; j++) {
				p += sizeof(int);
				p += sizeof(int);
				long sa_size = *(long*)(z + p); p += sizeof(long);
				int occ = sa_size < 500 ? sa_size : 500;
				p += occ * sizeof(long);
				zip_occ += occ;
			}
			free(z);
		}
		fprintf(stderr, "%ld %ld %ld %ld\n", mem_substring, zip_substring, mem_occ, zip_occ);
//		kt_for(16, worker1, data, data->n_seqs);
		free(data->zip_seeds); free(data->mem_seeds);
		aux->n_processed += data->n_seqs;
		double rtime_e = realtime();
//		fprintf(stderr, "[%s_step2] Check seeds for %d reads in %.1f real sec, have processed %ld reads in total\n",
//		        __func__, data->n_seqs, rtime_e-rtime_s, aux->n_processed);
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
	fprintf(stderr, "Seeds consistency: %.2f %% = %ld / %ld\n", 100.0 * p->hit_seed_n / p->mem_seed_n, p->hit_seed_n, p->mem_seed_n);
	fclose(aux.zip_seeds); fclose(aux.mem_seeds);
}

/************************
 * Check for alignment  *
 ************************/

#define SAMF_MUL_SEG    0x1   // template having multiple segments in sequencing
#define SAMF_BOTH_ALI   0x2   // each segment properly aligned according to the aligner
#define SAMF_UNMAP      0x4   // segment unmapped
#define SAMF_NEXT_UNMAP 0x8   // next segment in the template unmapped
#define SAMF_RC         0x10  // SEQ being reverse complemented
#define SAMF_NEXT_RC    0x20  // SEQ of the next segment in the template being reverse complemented
#define SAMF_READ1      0x40  // the first segment in the template
#define SAMF_READ2      0x80  // the last segment in the template
#define SAMF_SEC_ALI    0x100 // secondary alignment, multiple mapping
#define SAMF_N0_FLT     0x200 // not passing filters, such as platform/vendor quality controls
#define SAMF_PCR        0x400 // PCR or optical duplicate
#define SAMF_SUP_ALI    0x800 // supplementary alignment, chimeric alignment

static inline int mul_seg(int flag) { return ((flag & SAMF_MUL_SEG) != 0); }
static inline int both_ali(int flag) { return ((flag & SAMF_BOTH_ALI) != 0); }
static inline int unmap(int flag) { return ((flag & SAMF_UNMAP) != 0); }
static inline int next_unmap(int flag) { return ((flag & SAMF_NEXT_UNMAP) != 0); }
static inline int is_rc(int flag) { return ((flag & SAMF_RC) != 0); }
static inline int next_rc(int flag) { return ((flag & SAMF_NEXT_RC) != 0); }
static inline int is_read1(int flag) { return ((flag & SAMF_READ1) != 0); }
static inline int is_read2(int flag) { return ((flag & SAMF_READ2) != 0); }
static inline int sec_ali(int flag) { return ((flag & SAMF_SEC_ALI) != 0); }
static inline int no_flt(int flag) { return ((flag & SAMF_N0_FLT) != 0); }
static inline int pcr(int flag) { return ((flag & SAMF_PCR) != 0); }
static inline int sup_ali(int flag) { return ((flag & SAMF_SUP_ALI) != 0); }

static sam_line_t fetch_samline1(gzFile f) {
	char line[65536];
	sam_line_t sam; memset(&sam, 0, sizeof(sam));
	if (gzgets(f, line, sizeof(line)) == NULL) return sam;
	int len = (int)strlen(line);
	if(line[len-1] == '\n') { line[--len] = '\0'; }
	if (line[0] == '@') return sam; // Skip SAM header

	sam.data = strdup(line);
	int i, cnt_d = 0;
	for(i = 0; i < len; ++i) {
		if(sam.data[i] == '\t') {
			sam.data[i] = '\0';
		}
	}
	for(i = 0; i < len; ++i) {
		if(i == 0 || sam.data[i-1] == '\0') {
			if(cnt_d == 0) {
				sam.qname = sam.data + i;
			} else if(cnt_d == 1) {
				sam.flag = strtol(sam.data + i, NULL, 10);
			} else if(cnt_d == 2) {
				sam.rname = sam.data + i;
			} else if(cnt_d == 3) {
				sam.pos = strtol(sam.data + i, NULL, 10);
			} else if(cnt_d == 4) {
				sam.mapq = strtol(sam.data + i, NULL, 10);
			} else if(cnt_d == 5) {
				sam.cigar = sam.data + i;
			} else if(cnt_d == 6) {
				sam.rnext = sam.data + i;
			} else if(cnt_d == 7) {
				sam.pnext = strtol(sam.data + i, NULL, 10);
			} else if(cnt_d == 8) {
				sam.tlen = strtol(sam.data + i, NULL, 10);
			} else if(cnt_d == 9) {
				sam.seq = sam.data + i;
			} else if(cnt_d == 10) {
				sam.qual = sam.data + i;
			} else if (sam.data[i] == 'A' && sam.data[i+1] == 'S') {
				sam.as = strtol(sam.data + i + 5, NULL, 10);
			} else if (sam.data[i] == 'N' && sam.data[i+1] == 'M') {
				sam.nm = strtol(sam.data + i + 5, NULL, 10);
			}
		}
		if(sam.data[i] == '\0') {
			++cnt_d;
		}
	}
	return sam;
}

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
			if (m->mapq <= 3) continue;
			prof[0].mem_align_n++;
			if (!strcmp(z->rname, m->rname) && abs(z->pos - m->pos) <= 5) {
				prof[0].hit_p_n++;
			}
			if (z->nm == m->nm) prof[0].hit_nm_n++;
			if (abs(z->mapq - m->mapq) <= 5) prof[0].hit_mq_n++;
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

	int i;
	for (i = 1; i < 256; i++) prof[0].mem_align_n += prof[i].mem_align_n;
	for (i = 1; i < 256; i++) prof[0].hit_p_n += prof[i].hit_p_n;
	for (i = 1; i < 256; i++) prof[0].hit_nm_n += prof[i].hit_nm_n;
	for (i = 1; i < 256; i++) prof[0].hit_mq_n += prof[i].hit_mq_n;
	const prof_t *p = &prof[0];
	fprintf(stderr, "Primary: %ld\n", p->mem_align_n);
	fprintf(stderr, "POS:     %ld\n", p->hit_p_n);
	fprintf(stderr, "NM:      %ld\n", p->hit_nm_n);
	fprintf(stderr, "MAPQ     %ld\n", p->hit_mq_n);
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
