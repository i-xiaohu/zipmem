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

typedef kvec_t(int) int_v;

static inline int intv_len(bwtintv_t *p) {
	return (int)p->info - (p->info>>32);
}
typedef kvec_t(mem_seed_t) mem_seed_v;

cs_aux_t* cs_aux_init(int n, const bseq1_t *reads) {
	cs_aux_t *a;
	a = calloc(1, sizeof(cs_aux_t));
	int i, bases_n = 0;
	a->mem_loc = malloc(n * sizeof(int));
	for (i = 0; i < n; i++) {
		a->mem_loc[i] = bases_n;
		bases_n += reads[i].l_seq;
	}
	a->s = calloc(bases_n, sizeof(uint8_t));
	a->cnt[0] = calloc(bases_n, sizeof(int));
	a->cnt[1] = calloc(bases_n, sizeof(int));
	a->cnt[2] = calloc(bases_n, sizeof(int));
	a->cnt[3] = calloc(bases_n, sizeof(int));
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

typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	int n;
	bseq1_t *seqs;
	smem_aux_t **mem_aux;
	cs_aux_t *cs_aux;
	con_seq_v csv;
	int_v *misbuf;
} worker_t;

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
		cs.cnt[read->seq[i]][i]++;
	}
	return cs;
}

static int try_to_match_cs(con_seq_t *cs, int os, int len, const uint8_t *bases) {
	int i, ol = cs->n - os < len ?cs->n - os :len;
	int e = 0, m = (int)(0.05 * ol + 0.499);
	for (i = 0; i < ol; i++) {
		if (bases[i] != cs->s[os + i]) e++;
		if (e > m) return 0;
	}
	cs->n = cs->n > os + len ?cs->n :os + len;
	for (i = 0; i < len; i++) {
		if (bases[i] > 3) continue; // 'N' does not vote
		cs->cnt[bases[i]][os + i]++;
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
	w.chunks_n = threads_n; // Split the reads set into $threads_n uniform chunks
	w.n = n;
	w.reads = reads;
	w.aux = aux;
	w.chunk_csv = calloc(w.chunks_n, sizeof(con_seq_v));
	kt_for(threads_n, gen_cs_worker, &w, w.chunks_n);

	// Sum up CS of all chunks
	con_seq_v csv; kv_init(csv);
	int i, j;
	for(i = 0; i < w.chunks_n; i++) {
		int reads_l = (n+w.chunks_n-1)/w.chunks_n * i;
		int reads_r = reads_l + (n+w.chunks_n-1)/w.chunks_n;
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
static uint8_t* mem_sal(const mem_opt_t  *opt, const bwt_t *bwt, const bwtintv_v *mems) {
	int i;
	long bytes_n = sizeof(int);
	for (i = 0; i < mems->n; i++) { // Seed cluster
		bytes_n += sizeof(int) + sizeof(int) + sizeof(long); // qb, l, SA size,
		const bwtintv_t *p = &mems->a[i];
		bytes_n += sizeof(long) * (p->x[2] > opt->max_occ ? opt->max_occ : p->x[2]); // Occurrence locations
	}
	uint8_t *seeds = malloc(bytes_n);
	long pointer = 0;
	*(int*)(seeds + pointer) = mems->n; pointer += sizeof(int);
	for (i = 0; i < mems->n; ++i) {
		bwtintv_t *p = &mems->a[i];
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
	return seeds;
}

static void cs_seeding_worker(void *data, long cs_id, int t_id) {
	worker_t *w = (worker_t*)data;
	con_seq_t *cs = &w->csv.a[cs_id];
	const bwt_t *bwt = w->bwt;
	const mem_opt_t *opt = w->opt;
	smem_aux_t *aux = w->mem_aux[t_id];

	// BWA-MEM seeding on CS.
	mem_collect_intv(opt, bwt, cs->n, cs->s, aux);
	cs->seeds = mem_sal(opt, bwt, &aux->mem);
}

static void sup_seeding_worker(void *data, long seq_id, int t_id) {
	worker_t *w = (worker_t*)data;
	const bwt_t *bwt = w->bwt;

	const mem_opt_t *opt = w->opt;
	bseq1_t *read = &w->seqs[seq_id];
	int i, j, len = read->l_seq;
	uint8_t *bases = (uint8_t*)read->seq;
	int cs_id = w->cs_aux->belong_to[seq_id];
	assert(cs_id < w->csv.n);
	int cs_os = w->cs_aux->offset[seq_id];
	const con_seq_t *cs = &w->csv.a[cs_id];
	assert(cs_os + len <= cs->n);

	uint8_t is_rc = w->cs_aux->is_rc[seq_id];
	if (is_rc) reverse_complement(len, bases);

//	smem_aux_t *a = w->mem_aux[t_id];
	int_v *misbuf = &w->misbuf[t_id]; misbuf->n = 0;
	for (i = 0; i < len; i++) {
		if (bases[i] > 3) continue;
		if (cs->s[cs_os + i] != bases[i]) kv_push(int, *misbuf, i);
	}
	zmp.mis_n[t_id] += misbuf->n;

	// Reusing preparation
	long pointer = 0, bytes_n = 0;
	int csmem_n = *(int*)(cs->seeds + pointer); pointer += sizeof(int);
	int last_qb = -1;
	for (i = 0; i < csmem_n; i++) {
		int sb = *(int*)(cs->seeds + pointer); pointer += sizeof(int);
		int se = sb + *(int*)(cs->seeds + pointer); pointer += sizeof(int);
		long sa_size = *(long*)(cs->seeds + pointer); pointer += sizeof(long);
		pointer += ((sa_size < opt->max_occ) ?sa_size :opt->max_occ) * sizeof(long);
		assert(sb >= last_qb);
		last_qb = sb;
	}
	pointer = sizeof(int);

	for (i = 0; i < csmem_n; i++) {
		long mem_beg = pointer;
		int sb = *(int*)(cs->seeds + pointer); pointer += sizeof(int);
		if (sb >= cs_os + len) break;
		int se = sb + *(int*)(cs->seeds + pointer); pointer += sizeof(int);
		long sa_size = *(long*)(cs->seeds + pointer); pointer += sizeof(long);
		pointer += ((sa_size < opt->max_occ) ?sa_size :opt->max_occ) * sizeof(long);
		long mem_end = pointer;
		if (se <= cs_os) continue;
		bytes_n += mem_end - mem_beg;
	}

	long p2 = 0;
	int reused_mem = 0;
	read->sam = malloc(bytes_n * sizeof(uint8_t) + sizeof(long) + sizeof(int));
	p2 += sizeof(long); // Skip the reused mem bytes
	p2 += sizeof(int); // Skip the reused mem number
	if (bytes_n > 0) {
		// Reusing CS seeds
		pointer = sizeof(int);
		for (i = 0; i < csmem_n; i++) {
			int sb = *(int*)(cs->seeds + pointer); pointer += sizeof(int);
			if (sb >= cs_os + len) break;
			int se = sb + *(int*)(cs->seeds + pointer); pointer += sizeof(int);
			long sa_size = *(long*)(cs->seeds + pointer); pointer += sizeof(long);
			int true_occ = sa_size < opt->max_occ ?sa_size :opt->max_occ;
			// Drag the MEM down to read
			int l_clip = sb < cs_os ?cs_os - sb :0;
			int r_clip = se > cs_os+len ?se-cs_os-len :0;
			int down_sb = sb - cs_os + l_clip;
			int down_l = se - sb - l_clip - r_clip;
			if (se <= cs_os || down_l < opt->min_seed_len) {
				pointer += true_occ * sizeof(long);
				continue;
			}
			down_sb = is_rc ?len-down_sb-down_l :down_sb;
			reused_mem++;
			zmp.reused_n[t_id] += sa_size;
			*(int*)(read->sam + p2) = down_sb; p2 += sizeof(int);
			*(int*)(read->sam + p2) = down_l; p2 += sizeof(int);
			*(long*)(read->sam + p2) = sa_size; p2 += sizeof(long);
			for (j = 0; j < true_occ; j++) {
				long loc = *(long*)(cs->seeds + pointer); pointer += sizeof(long);
				loc += l_clip;
				loc = is_rc ?bwt->seq_len-loc-down_l :loc;
				*(long*)(read->sam + p2) = loc;  p2 += sizeof(long);
			}
		}
	}
	*(long*)read->sam = p2 - sizeof(long);
	*(int*)(read->sam + sizeof(long)) = reused_mem;



//	/* 1. Seeding */
//	int x, j, start_width = 1;
//	int can_reuse = 1; // can reuse then zip-seeding, otherwise throw it into BWA-MEM seeding.
//	misbuf->n = 0;
//	a->mem.n = 0;
//
//	const con_seq_t *b = &w->csv.a[cs_id];
//
//	/* Seeding */
//	if(!can_reuse) {
//		x = 0;
//		while (x < len) {
//			if (bases[x] < 4) {
//				x = bwt_smem1(bwt, len, bases, x, start_width, &a->mem1, a->tmpv);
//				for (i = 0; i < a->mem1.n; ++i) {
//					bwtintv_t *p = &a->mem1.a[i];
//					if (intv_len(p) >= opt->min_seed_len) {
//						kv_push(bwtintv_t, a->mem, *p);
//					}
//				}
//			} else ++x;
//		}
//	} else {
//		/* Reads in CS and not be kicked off, seeding on each mismatched position */
//		const con_seq_t *b = &w->csv.a[cs_id];
//		int last_x = -1;
//		for(i = 0; i < misbuf->n; ++i) {
//			int mis = misbuf->a[i];
//			if(mis < last_x) continue; // This mismatch has been covered by previous EMs.
//			if(bases[mis] > 3) continue; // Can't seeding on an ambiguous base N
//			int sum = b->cnt[0][dis+mis] + b->cnt[1][dis+mis] + b->cnt[2][dis+mis] + b->cnt[3][dis+mis];
//			if(b->cnt[bases[mis]][dis + mis] <= sum / 4) continue; // The mismatched base takes few percentage.
//			last_x = bwt_smem1(bwt, len, bases, mis, start_width, &a->mem1, a->tmpv);
//			for(j = 0; j < a->mem1.n; ++j) {
//				bwtintv_t *p = &a->mem1.a[j];
//				if(intv_len(p) >= opt->min_seed_len) {
//					kv_push(bwtintv_t, a->mem, *p);
//				}
//			}
//		}
//	}
//
//	/* 2. Re-seeding */
//	int old_n = (int)a->mem.n;
//	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
//	for(i = 0; i < old_n; ++i) {
//		bwtintv_t *p = &a->mem.a[i];
//		int start = (int)(p->info>>32), end = (int32_t)p->info;
//		if (end - start < split_len || p->x[2] > opt->split_width) continue;
//		bwt_smem1(bwt, len, bases, (start + end) >> 1, p->x[2] + 1, &a->mem1, a->tmpv);
//		for(j = 0; j < a->mem1.n; ++j) {
//			if(intv_len(&a->mem1.a[j]) >= opt->min_seed_len) {
//				kv_push(bwtintv_t, a->mem, a->mem1.a[j]);
//			}
//		}
//	}
//
//	/* 3. LAST-like */
//	if(can_reuse == 0 && opt->max_mem_intv > 0) {
//		x = 0;
//		while(x < len) {
//			if(bases[x] < 4) {
//				bwtintv_t m;
//				x = bwt_seed_strategy1(bwt, len, bases, x, opt->min_seed_len, opt->max_mem_intv, &m);
//				if(m.x[2] > 0) {
//					kv_push(bwtintv_t, a->mem, m);
//				}
//			} else {
//				++x;
//			}
//		}
//	}
//
//	if(can_reuse == 0) {
//		w->cs_id[seq_id] = -1; // Mark this for seeds sorting.
//		ks_introsort(mez_intv, a->mem.n, a->mem.a);
//	}
//
//	/* 4. SAL */
//	mem_seed_v *seeds = &w->seq_seeds[t_id]; seeds->n = 0;
//	for (i = 0; i < a->mem.n; ++i) {
//		bwtintv_t *p = &a->mem.a[i];
//		em2seeds(bwt, p, opt->max_occ, seeds);
//	}
//
//	if(can_reuse == 0) {
//		return ;
//	}
//
//	/* 5. Reuse shared seeds in CS */
//	const con_seq_t *b = &w->csv.a[cs_id];
//	const mem_seed_v *ss = &b->shared_seeds;
//	for(i = 0; i < ss->n; ++i) { // seeds with query-begin < dis + len
//		mem_seed_t *s = &ss->a[i];
//		if(s->qbeg >= dis + len || s->qbeg + s->len <= dis) { // has no intersection
//			continue;
//		}
//		int start = s->qbeg - dis, end = s->qbeg + s->len - dis; // [start, end) had fall down query read
//		int cut_l = (start > 0) ?start :0; // cut [start, end)
//		int cut_r = (end < len) ?end :len;
//		assert(cut_l < cut_r);
//		if(cut_r - cut_l < opt->min_seed_len) continue;
//		int_v *misv = &w->seedmis[t_id]; misv->n = 0;
//		kv_push(int, *misv, cut_l - 1);
//		for(j = 0; j < misbuf->n; ++j) {
//			int mis = misbuf->a[j];
//			if(mis >= cut_l && mis < cut_r) {
//				kv_push(int, *misv, mis);
//			}
//		}
//		kv_push(int, *misv, cut_r);
//		for(j = 0; j < misv->n - 1; ++j) {
//			int qb = misv->a[j] + 1, qe = misv->a[j+1]; // [qb, qe)
//			if(qe - qb < opt->min_seed_len) continue;
//			mem_seed_t rs;
//			rs.rbeg = s->rbeg + (qb - start);
//			rs.qbeg = qb;
//			rs.len = qe - qb;
//			rs.score = rs.len; // reused seeds didn't call em2seeds, so we need to assign score manually.
//			kv_push(mem_seed_t, *seeds, rs);
//		}
//	}
}

void zipmem_seeding(const mem_opt_t *opt, const bwt_t *bwt, int n, bseq1_t *seqs) {
	worker_t w;
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

	/* Seeding on Consensus Sequences */
	w.mem_aux = malloc(opt->n_threads * sizeof(smem_aux_t*));
	for (i = 0; i < opt->n_threads; i++) w.mem_aux[i] = smem_aux_init();
	kt_for(opt->n_threads, cs_seeding_worker, &w, w.csv.n);
	double csseed_c = cputime(), csseed_r = realtime();
	zmp.t_csseed[0] += csseed_c - gencs_c; zmp.t_csseed[1] += csseed_r - gencs_r;

	/* Reusing and supplementary seeding */
	w.misbuf = calloc(opt->n_threads, sizeof(int_v));
	kt_for(opt->n_threads, sup_seeding_worker, &w, n);
	for (i = 0; i < w.csv.n; i++) free(w.csv.a[i].seeds);
	cs_aux_destroy(w.cs_aux);
	free(w.csv.a);
	for (i = 0; i < opt->n_threads; i++) { smem_aux_destroy(w.mem_aux[i]); } free(w.mem_aux);
	for (i = 0; i < opt->n_threads; i++) { free(w.misbuf[i].a); } free(w.misbuf);
	double supseed_c = cputime(), supseed_r = realtime();
	zmp.t_supseed[0] += supseed_c - csseed_c; zmp.t_supseed[1] += supseed_r - csseed_r;
}
