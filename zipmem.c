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

typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	int n;
	bseq1_t *seqs;
	smem_aux_t **mem_aux;
	char **F1, **R1, **F2, **R2;
	int *offset;    // offset on CS of all reads
	uint8_t *is_rc; // if reads are reversed/complement in CS
	int *belong_to; // which CS the read belongs to
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

char complement[256];

static void gen_cs_worker(void *data, long seq_id, int t_id) {
	worker_t *w = (worker_t *)data;
	int *offset = w->offset;
	if (seq_id == 0) {
		offset[seq_id] = -1;
		return ;
	}
	const bseq1_t *r1 = &w->seqs[seq_id-1];
	const bseq1_t *r2 = &w->seqs[seq_id];
	char *F1 = w->F1[t_id], *R1 = w->R1[t_id];
	char *F2 = w->F2[t_id], *R2 = w->R2[t_id];
	uint8_t *is_rc = w->is_rc;
	int i, j, d;
	for (i = 0; i < r1->l_seq; i++) {
		F1[i] = r1->seq[i];
		R1[i] = complement[r1->seq[r1->l_seq-1-i]];
	}
	for (i = 0; i < r2->l_seq; i++) {
		F2[i] = r2->seq[i];
		R2[i] = complement[r2->seq[r2->l_seq-1-i]];
	}
	is_rc[seq_id] = 0;
	offset[seq_id] = -1;
	for (d = 0; d < r1->l_seq/2; d++) {
		int FF = 0, FR = 0, RF = 0, RR = 0;
		int ol = r1->l_seq - d < r2->l_seq ?r1->l_seq - d :r2->l_seq;
		int mis_bound = (int)(ol * 0.05 + 0.499);
		for (i = d, j = 0; i < r1->l_seq && j < r2->l_seq; i++, j++) {
			FF += F1[i] != F2[j];
			FR += F1[i] != R2[j];
			RF += R1[i] != F2[j];
			RR += R1[i] != R2[j];
			if (FF > mis_bound && FR > mis_bound && RF > mis_bound && RR > mis_bound) break;
		}
		if (FF <= mis_bound) is_rc[seq_id] = 1;
		else if (FR <= mis_bound) is_rc[seq_id] = 2;
		else if (RF <= mis_bound) is_rc[seq_id] = 3;
		else if (RR <= mis_bound) is_rc[seq_id] = 4;
		if (is_rc[seq_id]) {
			offset[seq_id] = d;
			break;
		}
	}
}

static con_seq_v connect_reads(int threads_n, worker_t *w) {
	int i;
	for (i = 0; i < 256; i++) complement[i] = (char)i;
	complement['A'] = 'T'; complement['a'] = 't';
	complement['C'] = 'G'; complement['c'] = 'g';
	complement['G'] = 'C'; complement['g'] = 'c';
	complement['T'] = 'A'; complement['t'] = 'a';

	int n = w->n, max_len = 0;
	const bseq1_t *reads = w->seqs;
	for (i = 0; i < n; i++) max_len = reads[i].l_seq > max_len ?reads[i].l_seq :max_len;
	w->F1 = malloc(sizeof(char*) * threads_n); w->R1 = malloc(sizeof(char*) * threads_n);
	w->F2 = malloc(sizeof(char*) * threads_n); w->R2 = malloc(sizeof(char*) * threads_n);
	for (i = 0; i < threads_n; i++) {
		w->F1[i] = malloc(max_len * sizeof(char)); w->R1[i] = malloc(max_len * sizeof(char));
		w->F2[i] = malloc(max_len * sizeof(char)); w->R2[i] = malloc(max_len * sizeof(char));
	}
	w->offset = malloc(n * sizeof(int));
	w->is_rc = malloc(n * sizeof(uint8_t));
	w->belong_to = malloc(n * sizeof(int));

	kt_for(threads_n, gen_cs_worker, w, n);

	for (i = 0; i < threads_n; i++) {
		free(w->F1[i]); free(w->R1[i]);
		free(w->F2[i]); free(w->R2[i]);
	}
	free(w->F1); free(w->R1);
	free(w->F2); free(w->R2);

	con_seq_v ret; kv_init(ret);
	con_seq_t tmp; tmp.l = 0;
	for (i = 1; i < n; i++) {
		if (w->offset[i] == -1) {
			tmp.r = i;
			kv_push(con_seq_t, ret, tmp);
			tmp.l = i;
		}
	}
	tmp.r = n;
	kv_push(con_seq_t, ret, tmp);
	zmp.reads_n += n;
	zmp.read_len = reads[0].l_seq;
	return ret;
}

// Look up SA for occurrence location on reference of MEMs.
static uint8_t* mem_sal(const mem_opt_t  *opt, const bwt_t *bwt, const bwtintv_v *mems) {
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
	return seeds;
}

static void cs_seeding(const mem_opt_t *opt, const bwt_t *bwt, con_seq_t *cs, const bseq1_t *reads,
					   int *offset, uint8_t *is_rc, smem_aux_t *aux) {
	int l = cs->l, r = cs->r;
	// Preparing
	int i, j, ml = 0, mr = reads[l].l_seq, p = 0;
	char dir = 'F'; // Assume the first read on F direction
	offset[l] = 0; is_rc[l] = 0;
	for (i = l+1; i < r; i++) {
		if (offset[i] == -1) break;
		if (is_rc[i] == 1) { // FF
			if (dir == 'F') {
				p = p + offset[i];
				dir = 'F';
			} else {
				p = p - offset[i];
				dir = 'R';
			}
		} else if (is_rc[i] == 2) { // FR
			if (dir == 'F') {
				p = p + offset[i];
				dir = 'R';
			} else {
				p = p - offset[i];
				dir = 'F';
			}
		} else if (is_rc[i] == 3) { // RF
			if (dir == 'R') {
				p = p + offset[i];
				dir = 'F';
			} else {
				p = p - offset[i];
				dir = 'R';
			}
		} else { // RR
			if (dir == 'R') {
				p = p + offset[i];
				dir = 'R';
			} else {
				p = p - offset[i];
				dir = 'F';
			}
		}
		offset[i] = p;
		is_rc[i] = dir == 'F' ?0 :1;
		ml = ml < p ?ml :p;
		mr = mr > p + reads[i].l_seq ?mr :p + reads[i].l_seq;
	}
	cs->n = mr - ml;
	for (i = l; i < r; i++) offset[i] -= ml;

	// Voting
	int16_t *cnt[4];
	for (i = 0; i < 4; i++) cnt[i] = calloc(cs->n, sizeof(int16_t));
	for (i = l; i < r; i++) {
		const bseq1_t *b = &reads[i];
		for (j = 0; j < b->l_seq; j++) b->seq[j] = (uint8_t)nst_nt4_table[b->seq[j]];
		if (is_rc[i]) reverse_complement(b->l_seq, b->seq);
		for (j = 0; j < b->l_seq; j++) {
			uint8_t c = b->seq[j];
			if (c < 4) cnt[c][offset[i]+j]++; // 'N' does not vote
		}
	}
	cs->s = malloc(cs->n * sizeof(uint8_t));
	for (i = 0; i < cs->n; i++) {
		int max_c = -1, max_id;
		for (j = 0; j < 4; j++) {
			if (cnt[j][i] > max_c) {
				max_c = cnt[j][i];
				max_id = j;
			}
		}
		cs->s[i] = max_id;
	}
	for (i = 0; i < 4; i++) free(cnt[i]);

	// BWA-MEM seeding on CS.
	mem_collect_intv(opt, bwt, cs->n, cs->s, aux);
	cs->seeds = mem_sal(opt, bwt, &aux->mem);
}

static inline int intv_len(const bwtintv_t *p) {
	return (int)p->info - (p->info>>32);
}

static void reuse_and_supplement(const mem_opt_t *opt, const bwt_t *bwt, const con_seq_t *cs,
								 bseq1_t *read, int cs_os, uint8_t is_rc, smem_aux_t *a, int t_id) {
	int i, j, len = read->l_seq;
	long total_occ = 0;

	// Reusing preparation
	long cs_p = 0, reused_bytes = 0;
	int csmem_n = *(int*)(cs->seeds + cs_p); cs_p += sizeof(int);
	for (i = 0; i < csmem_n; i++) {
		int sb = *(int*)(cs->seeds + cs_p); cs_p += sizeof(int);
		if (sb >= cs_os + len) break;
		int se = sb + *(int*)(cs->seeds + cs_p); cs_p += sizeof(int);
		long sa_size = *(long*)(cs->seeds + cs_p); cs_p += sizeof(long);
		int true_occ = sa_size < opt->max_occ ?sa_size :opt->max_occ;
		cs_p += true_occ * sizeof(long);
		int l_clip = sb < cs_os ?cs_os - sb :0;
		int r_clip = se > cs_os+len ?se-cs_os-len :0;
		int down_l = se - sb - l_clip - r_clip;
		if (se <= cs_os || down_l < opt->min_seed_len) continue;
		reused_bytes += sizeof(int) + sizeof(int) + sizeof(long) + true_occ * sizeof(long);
	}

	long reuse_p = 0;
	uint8_t *reused_mem = NULL;
	int reused_n = 0;
	if (reused_bytes > 0) {
		// Reusing CS seeds
		reused_mem = malloc(reused_bytes * sizeof(uint8_t));
		cs_p = sizeof(int);
		for (i = 0; i < csmem_n; i++) {
			int sb = *(int*)(cs->seeds + cs_p); cs_p += sizeof(int);
			if (sb >= cs_os + len) break;
			int se = sb + *(int*)(cs->seeds + cs_p); cs_p += sizeof(int);
			long sa_size = *(long*)(cs->seeds + cs_p); cs_p += sizeof(long);
			int true_occ = sa_size < opt->max_occ ?sa_size :opt->max_occ;
			// Drag the MEM down to read
			int l_clip = sb < cs_os ?cs_os - sb :0;
			int r_clip = se > cs_os+len ?se-cs_os-len :0;
			int down_sb = sb - cs_os + l_clip;
			int down_l = se - sb - l_clip - r_clip;
			if (se <= cs_os || down_l < opt->min_seed_len) {
				cs_p += true_occ * sizeof(long);
				continue;
			}
			down_sb = is_rc ?len-down_sb-down_l :down_sb;
			reused_n++;
			zmp.reused_n[t_id] += true_occ;
			*(int*)(reused_mem + reuse_p) = down_sb; reuse_p += sizeof(int);
			*(int*)(reused_mem + reuse_p) = down_l; reuse_p += sizeof(int);
			*(long*)(reused_mem + reuse_p) = sa_size; reuse_p += sizeof(long);
			total_occ += true_occ;
			for (j = 0; j < true_occ; j++) {
				long loc = *(long*)(cs->seeds + cs_p); cs_p += sizeof(long);
				loc += l_clip;
				loc = is_rc ?bwt->seq_len-loc-down_l :loc;
				*(long*)(reused_mem + reuse_p) = loc; reuse_p += sizeof(long);
			}
		}
	}
	assert(reuse_p == reused_bytes);

	long skip_bound = total_occ > 200 ?20 :200; // Memory protection
	uint8_t *bases = (uint8_t*)read->seq;

	// Supplementary seeding at each mismatched base
	a->mem.n = 0;
	int last_mem_end = -1;
	for (i = 0; i < len; i++) {
		if (bases[i] == cs->s[cs_os + i] || bases[i] > 3) continue;
		zmp.mis_n[t_id]++;
		if (i < last_mem_end) continue; // Skip the mismatch contained by previous mem
		last_mem_end = bwt_smem1(bwt, len, bases, i, 1, &a->mem1, a->tmpv);
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

	// Push the reused mems into sup mems array
	reuse_p = 0;
	for (i = 0; i < reused_n; i++) {
		bwtintv_t temp;
		temp.x[2] = 0; // Use x[2]=0 to mark reused mem
		temp.x[0] = reuse_p; // x[2] is used to record corresponding pointer
		int sb = *(int*)(reused_mem + reuse_p); reuse_p += sizeof(int);
		int se = sb + *(int*)(reused_mem + reuse_p); reuse_p += sizeof(int);
		long sa_size = *(long*)(reused_mem + reuse_p); reuse_p += sizeof(long);
		int true_occ = sa_size < opt->max_occ ?sa_size :opt->max_occ;
		reuse_p += true_occ * sizeof(long);
		temp.info = ((long)sb << 32) | se;
		kv_push(bwtintv_t, a->mem, temp);
	}

	// Sort the sup and reused mems together
	ks_introsort_mem_intv(a->mem.n, a->mem.a);

	long total_bytes = sizeof(long) + sizeof(int) + sup_bytes + reused_bytes;
	read->sam = malloc(total_bytes);
	long pointer = 0;
	pointer += sizeof(long); // Skip #bytes
	pointer += sizeof(int); // Skip #mems
	for (i = 0; i < a->mem.n; i++) {
		const bwtintv_t *p = &a->mem.a[i];
		if (p->x[2] == 0) {
			reuse_p = p->x[0];
			int sb = *(int*)(reused_mem + reuse_p); reuse_p += sizeof(int);
			int se = sb + *(int*)(reused_mem + reuse_p); reuse_p += sizeof(int);
			long sa_size = *(long*)(reused_mem + reuse_p);
			int true_occ = sa_size < opt->max_occ ?sa_size :opt->max_occ;
			long move_bytes = sizeof(int) + sizeof(int) + sizeof(long) + sizeof(long) * true_occ;
			memcpy(read->sam + pointer, reused_mem + p->x[0], move_bytes);
			pointer += move_bytes;
			assert(sb == p->info >> 32);
			assert(se == (int)p->info);
		} else {
			int count, qb = p->info>>32;
			int slen = (uint32_t)p->info - qb;
			qb = is_rc ?len-qb-slen :qb;
			*(int*)(read->sam + pointer) = qb; pointer += sizeof(int);
			*(int*)(read->sam + pointer) = slen; pointer += sizeof(int);
			*(long*)(read->sam + pointer) = p->x[2]; pointer += sizeof(long);
			int64_t k, step;
			step = p->x[2] > opt->max_occ? p->x[2] / opt->max_occ : 1;
			for (k = count = 0; k < p->x[2] && count < opt->max_occ; k += step, ++count) {
				int64_t rb = bwt_sa(bwt, p->x[0] + k);
				rb = is_rc ?bwt->seq_len-rb-slen :rb;
				*(long*)(read->sam + pointer) = rb; pointer += sizeof(long);
			}
		}
	}
	*(long*)read->sam = total_bytes - sizeof(long);
	*(int*)(read->sam + sizeof(long)) = a->mem.n;
	free(reused_mem);
}

static void seeding_worker(void *data, long cs_id, int t_id) {
	worker_t *w = (worker_t*)data;
	con_seq_t *cs = &w->csv.a[cs_id];
	const bwt_t *bwt = w->bwt;
	const mem_opt_t *opt = w->opt;
	smem_aux_t *aux = w->mem_aux[t_id];
	bseq1_t *reads = w->seqs;

	cs_seeding(opt, bwt, cs, reads, w->offset, w->is_rc, aux);

	int i, j;
	zmp.cs_bases += cs->n;
	for (i = cs->l; i < cs->r; i++) {
		const bseq1_t *b = &reads[i];
		for (j = 0; j < b->l_seq; j++) {
			if (b->seq[j] != cs->s[w->offset[i] + j]) zmp.mis_n[t_id]++;
		}
	}

	for (i = cs->l; i < cs->r; i++) {
		reuse_and_supplement(opt, bwt, cs, &reads[i], w->offset[i], w->is_rc[i], aux, t_id);
	}
	free(cs->seeds); free(cs->s);
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
	w.csv = connect_reads(opt->n_threads, &w);
	double gencs_c = cputime(), gencs_r = realtime();
	zmp.t_gencs[0] += gencs_c - ctime; zmp.t_gencs[1] += gencs_r - rtime;

	/* Seeding on Consensus Sequences */
	/* Reusing and supplementary seeding */
	w.mem_aux = malloc(opt->n_threads * sizeof(smem_aux_t*));
	for (i = 0; i < opt->n_threads; i++) w.mem_aux[i] = smem_aux_init();
	kt_for(opt->n_threads, seeding_worker, &w, w.csv.n);
	for (i = 0; i < opt->n_threads; i++) { smem_aux_destroy(w.mem_aux[i]); } free(w.mem_aux);
	free(w.offset); free(w.is_rc); free(w.csv.a);
	double csseed_c = cputime(), csseed_r = realtime();
	zmp.t_csseed[0] += csseed_c - gencs_c; zmp.t_csseed[1] += csseed_r - gencs_r;
}
