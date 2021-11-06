//
// Created by ixiaohu on 2021/9/21.
//

#include <limits.h>

#include "cstl//kvec.h"
#include "cstl//ksort.h"
#include "cstl/kbtree.h"
#include "cstl/kthread.h"
#include "bwalib/utils.h"
#include "time_prof.h"
#include "zipmem.h"

extern zsmem_prof_t zmp;

cs_aux_t* cs_aux_init(int n, const bseq1_t *reads) {
	// For 160M bases, ~820M bytes needed.
	cs_aux_t *a;
	a = calloc(1, sizeof(cs_aux_t));
	int i, bases_n = 0;
	a->mem_loc = malloc(n * sizeof(int));
	for (i = 0; i < n; i++) {
		a->mem_loc[i] = bases_n;
		bases_n += reads[i].l_seq;
	}
	a->s = calloc(bases_n, sizeof(uint8_t));
	a->cnt[0] = calloc(bases_n, sizeof(uint8_t));
	a->cnt[1] = calloc(bases_n, sizeof(uint8_t));
	a->cnt[2] = calloc(bases_n, sizeof(uint8_t));
	a->cnt[3] = calloc(bases_n, sizeof(uint8_t));
	a->offset = malloc(n * sizeof(int));
	a->is_rc = malloc(n * sizeof(uint8_t));
	a->belong_to = malloc(n * sizeof(int));
	return a;
}

void cs_aux_destroy(cs_aux_t *a)  {
	free(a->mem_loc);
	free(a->s);
	free(a->cnt[0]);
	free(a->cnt[1]);
	free(a->cnt[2]);
	free(a->cnt[3]);
	free(a->offset);
	free(a->is_rc);
	free(a->belong_to);
	free(a);
}

static inline void reverse_complement(int len, uint8_t *bases) {
	int i;
	for(i = 0; i < len; i++) {
		bases[i] = bases[i] > 3 ? bases[i] : (3 - bases[i]);
	}
	for (i = 0; i < len/2; i++) {
		uint8_t temp = bases[i];
		bases[i] = bases[len-1-i];
		bases[len-1-i] = temp;
	}
}

typedef struct {
	int chunks_n;
	int n;
	bseq1_t *reads;
	cs_aux_t *aux;
	con_seq_v *chunk_csv;
} cs_worker_t;

static con_seq_t init_cs(int read_id, const bseq1_t *read, cs_aux_t *aux) {
	con_seq_t cs;
	cs.s = aux->s + aux->mem_loc[read_id];
	cs.cnt[0] = aux->cnt[0] + aux->mem_loc[read_id];
	cs.cnt[1] = aux->cnt[1] + aux->mem_loc[read_id];
	cs.cnt[2] = aux->cnt[2] + aux->mem_loc[read_id];
	cs.cnt[3] = aux->cnt[3] + aux->mem_loc[read_id];
	cs.n = read->l_seq;
	int i;
	for (i = 0; i < read->l_seq; i++) {
		if (read->seq[i] > 3) continue;
		cs.s[i] = read->seq[i];
		if (cs.cnt[read->seq[i]][i] < UINT8_MAX) cs.cnt[read->seq[i]][i]++;
	}
	return cs;
}

static inline int try_to_match_cs(con_seq_t *cs, int os, int len, const uint8_t *bases) {
	int i, ol = cs->n - os < len ?cs->n - os :len;
	int e = 0, m = (int)(0.05 * ol + 0.499);
	for (i = 0; i < ol; i++) {
		if (bases[i] != cs->s[os + i]) e++;
		if (e > m) return 0;
	}
	cs->n = cs->n > os + len ?cs->n :os + len;
	for (i = 0; i < len; i++) {
		if (bases[i] > 3) continue; // 'N' does not vote
		if (cs->cnt[bases[i]][os + i] < UINT8_MAX) cs->cnt[bases[i]][os + i]++;
		else continue;
		if (cs->s[os + i] == bases[i]) continue;
		if (cs->cnt[bases[i]][os + i] > cs->cnt[cs->s[os + i]][os + i]) {
			cs->s[os + i] = bases[i];
		}
	}
	return 1;
}

static void gen_cs_worker(void *data, long c_id, int t_id) {
	cs_worker_t *w = (cs_worker_t *)data;
	bseq1_t *reads = w->reads;
	cs_aux_t *aux = w->aux;

	con_seq_v *csv = &w->chunk_csv[c_id]; kv_init(*csv);
	int *offset = w->aux->offset;
	uint8_t *is_rc = w->aux->is_rc;
	int *belong_to = w->aux->belong_to;
	int reads_l = (w->n+w->chunks_n-1)/w->chunks_n * c_id;
	int reads_r = reads_l + (w->n+w->chunks_n-1)/w->chunks_n;
	reads_r = reads_r < w->n ?reads_r :w->n;

	int i, j, d;
	for (i = reads_l; i < reads_r; i++) {
		int l_seq = reads[i].l_seq;
		uint8_t *seq = (uint8_t*)reads[i].seq;
		for(j = 0; j < l_seq; ++j) {
			seq[j] = seq[j] < 5 ?seq[j] :nst_nt4_table[seq[j]];
		}
	}
	// Take the first read to initialize a CS
	con_seq_t tcs = init_cs(reads_l, &reads[reads_l], aux);
	offset[reads_l] = 0;
	is_rc[reads_l] = 0;
	belong_to[reads_l] = csv->n;
	for(i = reads_l + 1; i < reads_r; i++) {
		int l_seq = reads[i].l_seq;
		uint8_t *seq = (uint8_t*)reads[i].seq;
		offset[i] = -1; is_rc[i] = 0;
		// enumerate matched position on CS
		for(d = offset[i-1]; d < tcs.n; d++) {
			if (tcs.n - d < l_seq / 2) break; // At least overlapping a half of read length
			int added = try_to_match_cs(&tcs, d, l_seq, seq);
			if (added) {
				offset[i] = d;
				belong_to[i] = csv->n;
				break;
			}
		}
		if (offset[i] != -1) continue;
		// Forward connecting failed, try RC
		is_rc[i] = 1;
		reverse_complement(l_seq, seq);
		for(d = offset[i-1]; d < tcs.n; d++) {
			if (tcs.n - d < l_seq / 2) break; // At least overlapping a half of read length
			int added = try_to_match_cs(&tcs, d, l_seq, seq);
			if (added) {
				offset[i] = d;
				belong_to[i] = csv->n;
				break;
			}
		}
		reverse_complement(l_seq, seq);
		if (offset[i] != -1) continue;
		// Both forward and RC connecting failed, start a new CS by $reads[i]
		kv_push(con_seq_t, *csv, tcs);
		tcs = init_cs(i, &reads[i], aux);
		offset[i] = 0;
		is_rc[i] = 0;
		belong_to[i] = csv->n;
	}
	kv_push(con_seq_t, *csv, tcs);
}

con_seq_v connect_reads_to_cs(int threads_n, int n, bseq1_t *reads, cs_aux_t *aux) {
	cs_worker_t w;
	int each_chunk_size = (n+threads_n-1) / threads_n; // ceil(n / threads)
	w.chunks_n = (n+each_chunk_size-1)/each_chunk_size; // Split the reads set into $threads_n uniform chunks
	w.n = n;
	w.reads = reads;
	w.aux = aux;
	w.chunk_csv = calloc(w.chunks_n, sizeof(con_seq_v));
	kt_for(threads_n, gen_cs_worker, &w, w.chunks_n);

	// Sum up CS of all chunks
	con_seq_v csv; kv_init(csv);
	int i, j;
	for(i = 0; i < w.chunks_n; i++) {
		int reads_l = each_chunk_size * i;
		int reads_r = reads_l + each_chunk_size;
		reads_r = reads_r < n ?reads_r :n;
		for(j = reads_l; j < reads_r; j++) {
			w.aux->belong_to[j] += csv.n;
		}
		const con_seq_v *p = &w.chunk_csv[i];
		for(j = 0; j < p->n; j++) {
			kv_push(con_seq_t, csv, p->a[j]);
		}
		free(p->a);
	}

	do {
		zmp.reads_n += n;
		zmp.read_len = reads[0].l_seq;
		for (i = 0; i < csv.n; i++) {
			con_seq_t *b = &csv.a[i];
			zmp.cs_bases += b->n;
		}
	} while (0);

	free(w.chunk_csv);
	return csv;
}

// Look up SA for occurrence location on reference of MEMs.
size_t cs_vram[256], read_vram[256], sup_vram[256];
static uint8_t* mem_sal(const mem_opt_t  *opt, const bwt_t *bwt, const bwtintv_v *mems, int t_id) {
	int i, cnt = 0;
	long bytes_n = sizeof(int);
	int occ_bound = 50; // We further restrict the repetitive level of CS seeds
                         // Because reusing repetitive seeds leads to large memory consumption.
	for (i = 0; i < mems->n; i++) { // Seed cluster
		const bwtintv_t *p = &mems->a[i];
		if (p->x[2] > occ_bound) continue; // Ignore the seeds with too many occurrences
		cnt++;
		bytes_n += sizeof(int) + sizeof(int) + sizeof(long); // qb, l, SA size,
		bytes_n += sizeof(long) * (p->x[2] > opt->max_occ ? opt->max_occ : p->x[2]); // Occurrence locations
	}
	uint8_t *seeds = malloc(bytes_n);
	long pointer = 0;
	*(int*)(seeds + pointer) = cnt; pointer += sizeof(int);
	for (i = 0; i < mems->n; ++i) {
		bwtintv_t *p = &mems->a[i];
		if (p->x[2] > occ_bound) continue;
		int count, slen = (uint32_t)p->info - (p->info>>32); // seed length
		*(int*)(seeds + pointer) = p->info>>32; pointer += sizeof(int);
		*(int*)(seeds + pointer) = slen; pointer += sizeof(int);
		*(long*)(seeds + pointer) = p->x[2]; pointer += sizeof(long);
		int64_t k, step;
		step = p->x[2] > opt->max_occ? p->x[2] / opt->max_occ : 1;
		for (k = count = 0; k < p->x[2] && count < opt->max_occ; k += step, ++count) {
			int64_t rb = bwt_sa(bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
			*(long*)(seeds + pointer) = rb; pointer += sizeof(long);
		}
	}
	assert(pointer == bytes_n);
	cs_vram[t_id] += bytes_n;
	return seeds;
}

static void cs_seeding_worker(void *data, long cs_id, int t_id) {
	zip_worker_t *w = (zip_worker_t *)data;
	con_seq_t *cs = &w->csv.a[cs_id];
	const bwt_t *bwt = w->bwt;
	const mem_opt_t *opt = w->opt;
	smem_aux_t *aux = w->mem_aux[t_id];

	// BWA-MEM seeding on CS.
	mem_collect_intv(opt, bwt, cs->n, cs->s, aux);
	cs->seeds = mem_sal(opt, bwt, &aux->mem, t_id);
}

static inline int intv_len(const bwtintv_t *p) {
	return (int)p->info - (p->info>>32);
}

static void sup_seeding_worker(void *data, long seq_id, int t_id) {
	zip_worker_t *w = (zip_worker_t *)data;
	const bwt_t *bwt = w->bwt;

	const mem_opt_t *opt = w->opt;
	bseq1_t *read = &w->seqs[seq_id];
	int i, j, len = read->l_seq;
	int cs_id = w->cs_aux->belong_to[seq_id];
	int cs_os = w->cs_aux->offset[seq_id];
	const con_seq_t *cs = &w->csv.a[cs_id];
	uint8_t is_rc = w->cs_aux->is_rc[seq_id];

	// Calculating the reused memory
	long pointer = 0, reuse_bytes = 0;
	int csmem_n = *(int*)(cs->seeds + pointer); pointer += sizeof(int);
	int reuse_n = 0;
	long total_occ = 0;
	for (i = 0; i < csmem_n; i++) {
		int sb = *(int*)(cs->seeds + pointer); pointer += sizeof(int);
		if (sb >= cs_os + len) break;
		int se = sb + *(int*)(cs->seeds + pointer); pointer += sizeof(int);
		long sa_size = *(long*)(cs->seeds + pointer); pointer += sizeof(long);
		int true_occ = sa_size < opt->max_occ ?sa_size :opt->max_occ;
		pointer += true_occ * sizeof(long);
		int l_clip = sb < cs_os ?cs_os - sb :0;
		int r_clip = se > cs_os+len ?se-cs_os-len :0;
		int down_l = se - sb - l_clip - r_clip;
		if (se <= cs_os || down_l < opt->min_seed_len) continue;
		reuse_n++;
		total_occ += true_occ;
		reuse_bytes += sizeof(int) + sizeof(int) + sizeof(long) + true_occ * sizeof(long);
	}
	zmp.reused_n[t_id] += total_occ;

	// Supplementary seeding
	long skip_bound = total_occ > 200 ?20 :200; // Memory protection
	uint8_t *bases = (uint8_t*)read->seq;
	if (is_rc) reverse_complement(len, bases);

	smem_aux_t *a = w->mem_aux[t_id]; a->mem.n = 0;
	int_v *misbuf = &w->misbuf[t_id]; misbuf->n = 0;
	for (i = 0; i < len; i++) {
		if (bases[i] > 3) continue;
		if (cs->s[cs_os + i] != bases[i]) {
			kv_push(int, *misbuf, i);
		}
	}
	zmp.mis_n[t_id] += misbuf->n;
	int last_mem_end = -1;
	for (i = 0; i < misbuf->n; i++) {
		int x = misbuf->a[i];
		if (x < last_mem_end) continue; // Skip the mismatch contained by previous mem
		last_mem_end = bwt_smem1(bwt, len, bases, x, 1, &a->mem1, a->tmpv);
		for (j = 0; j < a->mem1.n; j++) {
			const bwtintv_t *p = &a->mem1.a[j];
			if (p->x[2] > skip_bound) continue;
			if (intv_len(p) >= opt->min_seed_len) {
				kv_push(bwtintv_t, a->mem, *p);
			}
		}
	}

	/* 2. Re-seeding */
	int old_n = (int)a->mem.n;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	for(i = 0; i < old_n; i++) {
		const bwtintv_t *p = &a->mem.a[i];
		int start = (int)(p->info>>32), end = (int32_t)p->info;
		if (end - start < split_len || p->x[2] > opt->split_width) continue;
		bwt_smem1(bwt, len, bases, (start + end) / 2, p->x[2] + 1, &a->mem1, a->tmpv);
		for(j = 0; j < a->mem1.n; ++j) {
			if (a->mem1.a[j].x[2] > skip_bound) continue;
			if(intv_len(&a->mem1.a[j]) >= opt->min_seed_len) {
				kv_push(bwtintv_t, a->mem, a->mem1.a[j]);
			}
		}
	}

	// Preparing memory for supplementary seeds
	long sup_bytes = 0;
	for (i = 0; i < a->mem.n; i++) {
		const bwtintv_t *p = &a->mem.a[i];
		int true_occ = p->x[2] < opt->max_occ ?p->x[2] :opt->max_occ;
		zmp.sup_n[t_id] += true_occ;
		sup_bytes += sizeof(int) + sizeof(int) + sizeof(long) + true_occ * sizeof(long);
	}
	sup_vram[t_id] += sup_bytes;

	read->sam = malloc(sizeof(long) + sizeof(int) + sup_bytes);
	long p3 = 0;
	p3 += sizeof(long); // Skip #bytes
	p3 += sizeof(int); // Skip #mems
	for (i = 0; i < a->mem.n; i++) {
		const bwtintv_t *p = &a->mem.a[i];
		int count, qb = p->info>>32;
		int slen = (uint32_t)p->info - qb;
		qb = is_rc ?len-qb-slen :qb;
		*(int*)(read->sam + p3) = qb; p3 += sizeof(int);
		*(int*)(read->sam + p3) = slen; p3 += sizeof(int);
		*(long*)(read->sam + p3) = p->x[2]; p3 += sizeof(long);
		int64_t k, step;
		step = p->x[2] > opt->max_occ? p->x[2] / opt->max_occ : 1;
		for (k = count = 0; k < p->x[2] && count < opt->max_occ; k += step, ++count) {
			int64_t rb = bwt_sa(bwt, p->x[0] + k);
			rb = is_rc ?bwt->seq_len-rb-slen :rb;
			*(long*)(read->sam + p3) = rb; p3 += sizeof(long);
		}
	}
	*(long*)read->sam = sizeof(int) + reuse_bytes + sup_bytes;
	*(int*)(read->sam + sizeof(long)) = a->mem.n + reuse_n;
	assert(p3 == sizeof(long) + sizeof(int) + sup_bytes);
	if (is_rc) reverse_complement(len, bases);
}

zip_worker_t *zipmem_seeding(const mem_opt_t *opt, const bwt_t *bwt, int n, bseq1_t *seqs) {
	zip_worker_t w;
	double ctime, rtime;
	int i;

	ctime = cputime(); rtime = realtime();
	w.opt = opt;
	w.bwt = bwt;
	w.n = n;
	w.seqs = seqs;

	/* Get consensus sequences, parallel process each chunk */
	w.cs_aux = cs_aux_init(n, seqs);
	w.csv = connect_reads_to_cs(opt->n_threads, n, seqs, w.cs_aux);
	double gencs_c = cputime(), gencs_r = realtime();
	zmp.t_gencs[0] += gencs_c - ctime; zmp.t_gencs[1] += gencs_r - rtime;

	memset(cs_vram, 0, sizeof(cs_vram));
	memset(read_vram, 0, sizeof(read_vram));
	memset(sup_vram, 0, sizeof(sup_vram));

	/* Seeding on Consensus Sequences */
	w.mem_aux = malloc(opt->n_threads * sizeof(smem_aux_t*));
	for (i = 0; i < opt->n_threads; i++) w.mem_aux[i] = smem_aux_init();
	kt_for(opt->n_threads, cs_seeding_worker, &w, w.csv.n);
	double csseed_c = cputime(), csseed_r = realtime();
	zmp.t_csseed[0] += csseed_c - gencs_c; zmp.t_csseed[1] += csseed_r - gencs_r;
	for (i = 0; i < opt->n_threads; i++) { smem_aux_destroy(w.mem_aux[i]); } free(w.mem_aux);

	/* Reusing and supplementary seeding */
	w.mem_aux = malloc(opt->n_threads * sizeof(smem_aux_t*));
	for (i = 0; i < opt->n_threads; i++) w.mem_aux[i] = smem_aux_init();
	w.misbuf = calloc(opt->n_threads, sizeof(int_v));
	kt_for(opt->n_threads, sup_seeding_worker, &w, n);
//	for (i = 0; i < w.csv.n; i++) free(w.csv.a[i].seeds);
//	cs_aux_destroy(w.cs_aux);
//	free(w.csv.a);
	for (i = 0; i < opt->n_threads; i++) { smem_aux_destroy(w.mem_aux[i]); } free(w.mem_aux);
	for (i = 0; i < opt->n_threads; i++) { free(w.misbuf[i].a); } free(w.misbuf);
	double supseed_c = cputime(), supseed_r = realtime();
	zmp.t_supseed[0] += supseed_c - csseed_c; zmp.t_supseed[1] += supseed_r - csseed_r;

	for(i = 1; i < 256; i++) cs_vram[0] += cs_vram[i];
	for(i = 1; i < 256; i++) read_vram[0] += read_vram[i];
	for(i = 1; i < 256; i++) sup_vram[0] += sup_vram[i];

	zmp.cs_vram = zmp.cs_vram > cs_vram[0] ?zmp.cs_vram :cs_vram[0];
	zmp.read_vram = zmp.read_vram > read_vram[0] ?zmp.read_vram :read_vram[0];
	zmp.sup_vram = zmp.sup_vram > sup_vram[0] ?zmp.sup_vram :sup_vram[0];

	zip_worker_t *ret = malloc(sizeof(zip_worker_t));
	*ret = w;
	return ret;
}
