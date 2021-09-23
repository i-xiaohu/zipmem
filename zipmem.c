//
// Created by ixiaohu on 2021/9/21.
//

#include <limits.h>

#include "cstl//kvec.h"
#include "cstl//ksort.h"
#include "cstl/kbtree.h"
#include "bwalib/utils.h"
#include "time_prof.h"
#include "zipmem.h"
extern zsmem_prof_t zmp;

typedef kvec_t(int) int_v;

typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

static smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	a = calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = calloc(1, sizeof(bwtintv_v));
	return a;
}

static void smem_aux_destroy(smem_aux_t *a)
{
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

static inline int intv_len(bwtintv_t *p) {
	return (int)p->info - (p->info>>32);
}
typedef kvec_t(mem_seed_t) mem_seed_v;

typedef struct {
	int l, r; // [l, r)
	int n; // consensus sequence length
	uint8_t *s; // consensus sequence generated form the reads batch
	mem_seed_v shared_seeds; // seeds generated from s
	int *cnt[4]; // Saving votes.
} r_batch_t;
typedef kvec_t(r_batch_t) r_batch_v;

typedef struct {
	int m;
	int *cnt[4]; // counts of A,C,G,T in a reads batch
} cs_aux_t;

static cs_aux_t* cs_aux_init() {
	cs_aux_t *a;
	a = calloc(1, sizeof(cs_aux_t));
	a->m = 65536; // 64K * 4 * 4 = 1M bytes each thread.
	a->cnt[0] = malloc(a->m * sizeof(int));
	a->cnt[1] = malloc(a->m * sizeof(int));
	a->cnt[2] = malloc(a->m * sizeof(int));
	a->cnt[3] = malloc(a->m * sizeof(int));
	return a;
}

static void cs_aux_destroy(cs_aux_t *a) {
	free(a->cnt[0]);
	free(a->cnt[1]);
	free(a->cnt[2]);
	free(a->cnt[3]);
	free(a);
}

typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	smem_aux_t **aux;
	cs_aux_t **cs_aux;
	mem_seed_v *seq_seeds; // seeds temporary buffer.
	int_v *misbuf;
	int_v *seedmis;
	bseq1_t *seqs;
	int_v chunk_intv; // Chunk interval.
	uint8_t *is_rev; // The read maybe reversed in CS generating stage.
	r_batch_v *back_batches;
	r_batch_v batches;
	int n;
	int *dis;
	int *cs_id;
	int64_t n_processed;
} worker_t;

static inline void reverse(char *des, const char *src, int len) {
	int i;
	for(i = 0; i < len; ++i) {
		des[len-1-i] = src[i] > 3 ?src[i] :(3 - src[i]);
	}
}

static void gen_cs(void *data, int c_id, int t_id) {
	worker_t *w = (worker_t*)data;
	const mem_opt_t *opt = w->opt;
	bseq1_t *seqs = w->seqs;
	r_batch_v *batches = &w->back_batches[t_id];
	int *dis = w->dis;
	int reads_l = w->chunk_intv.a[c_id];
	int reads_r = w->chunk_intv.a[c_id + 1];
	int i, j, d, k, max_errors = 8, min_matches = 20, ref_l = seqs[reads_l].l_seq;
	char *ref = malloc(ref_l * sizeof(char));
	char *rev = malloc(ref_l * sizeof(char));
	int *votes[4];
	for(i = 0; i < 4; ++i) {
		votes[i] = malloc(ref_l * sizeof(int));
	}
	r_batch_t tmp; tmp.l = reads_l;
	for(i = reads_l; i < reads_r; ++i) {
		char *seq = seqs[i].seq;
		int l_seq = seqs[i].l_seq;
		dis[i] = INT_MAX;
		for(j = 0; j < l_seq; ++j) {
			seq[j] = nst_nt4_table[(int)seq[j]];
		}
		if(i > reads_l) {
			// try to connect to ref
			for(d = 0; d < ref_l - min_matches; ++d) { // enumerate distance
				int e = 0, m = 0;
				for(k = 0; k < l_seq && d + k < ref_l; ++k) { // enumerate bases
					// N is considered as a mismatch
					if(seq[k] != ref[d + k]) {
						++e;
						if(e > max_errors) {
							break;
						}
					} else {
						++m;
					}
				}
				if(e <= max_errors && m >= min_matches) {
					dis[i] = d;
					break;
				}
			}
			if(dis[i] == INT_MAX) {
				reverse(rev, seq, l_seq);
				memcpy(seq, rev, l_seq * sizeof(char));
				// try to connect to ref
				for(d = 0; d < ref_l - min_matches; ++d) { // enumerate distance
					int e = 0, m = 0;
					for(k = 0; k < l_seq && d + k < ref_l; ++k) { // enumerate bases
						if(seq[k] != ref[d + k]) {
							++e;
							if(e > max_errors) {
								break;
							}
						} else {
							++m;
						}
					}
					if(e <= max_errors && m >= min_matches) {
						w->is_rev[i] = 1;
						dis[i] = d;
						break;
					}
				}
				if(dis[i] == INT_MAX) {
					reverse(rev, seq, l_seq);
					memcpy(seq, rev, l_seq * sizeof(char));
				}
			}
		}
		if(dis[i] == INT_MAX) {
			// reset ref
			for(j = 0; j < 4; ++j) {
				memset(votes[j], 0, ref_l * sizeof(int));
			}
			if(i > reads_l) {
				tmp.r = i;
				if(tmp.r - tmp.l >= 2) {
					kv_push(r_batch_t, *batches, tmp);
				}
				tmp.l = i;
			}
		}
		// delete dis[i] bases front
		if(dis[i] != INT_MAX) {
			for(j = 0; j < ref_l - dis[i]; ++j) {
				ref[j] = ref[j + dis[i]];
				for(k = 0; k < 4; ++k) {
					votes[k][j] = votes[k][j + dis[i]];
				}
			}
			for(j = ref_l - dis[i]; j < ref_l; ++j) {
				for(k = 0; k < 4; ++k) {
					votes[k][j] = 0;
				}
			}
		}
		// voting of read(i), update the ref.
		for(j = 0; j < ref_l; ++j) {
			if(seq[j] > 3) continue; // N doesn't vote.
			++votes[(int)seq[j]][j];
		}
		for(j = 0; j < ref_l; ++j) {
			int max = -1, max_id = 0;
			for(k = 0; k < 4; ++k) {
				if(votes[k][j] > max) {
					max = votes[k][j];
					max_id = k;
				}
			}
			ref[j] = (uint8_t)max_id;
		}
	}
	tmp.r = reads_r;
	if(tmp.r - tmp.l >= 2) {
		kv_push(r_batch_t, *batches, tmp);
	}
	free(ref); free(rev);
	for(i = 0; i < 4; ++i) {
		free(votes[i]);
	}
}

static void em2seeds(const bwt_t *bwt, bwtintv_t *p, int max_occ, mem_seed_v *seeds) {
	int step, count, slen = intv_len(p);
	int64_t k;
	step = p->x[2] > max_occ? p->x[2] / max_occ : 1;
	for(k = count = 0; k < p->x[2] && count < max_occ; k += step, ++count) {
		mem_seed_t tmp;
		tmp.rbeg = bwt_sa(bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
		tmp.len = slen;
		tmp.score = slen;
		tmp.qbeg = (int)(p->info>>32);
		kv_push(mem_seed_t, *seeds, tmp);
	}
}

// Sort seeds by [qb, qe, rb]
#define mem_seed_lt_step2(a, b) ((a).len != (b).len ?(a).len < (b).len :(a).rbeg < (b).rbeg)
#define mem_seed_lt(a, b) ((a).qbeg != (b).qbeg ?(a).qbeg < (b).qbeg :mem_seed_lt_step2((a), (b)))
KSORT_INIT(mem_seeds, mem_seed_t, mem_seed_lt)

// Sort EMs by [qb, qe]
#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mez_intv, bwtintv_t, intv_lt)

static void cs_seeding(void *data, int cs_id, int t_id) {

	worker_t *w = (worker_t*)data;
	r_batch_t *b = &w->batches.a[cs_id];
	cs_aux_t *cs = w->cs_aux[t_id];
	const bwt_t *bwt = w->bwt;
	const mem_opt_t *opt = w->opt;
	smem_aux_t *a = w->aux[t_id];

	int i, j;
	b->n = 0;
	for(i = b->l; i < b->r; ++i) {
		w->cs_id[i] = cs_id;
		w->dis[i] = (i == b->l) ?0 :(w->dis[i-1] + w->dis[i]);
		b->n = (w->dis[i] + w->seqs[i].l_seq > b->n) ?(w->dis[i] + w->seqs[i].l_seq) :b->n;
	}

	/* Majority rule: voting for Consensus Sequence */
	b->s = malloc(b->n * sizeof(uint8_t));
	if(b->n > cs->m) {
		cs->m = b->n << 1;
		for(i = 0; i < 4; ++i) {
			cs->cnt[i] = realloc(cs->cnt[i], cs->m * sizeof(int));
		}
	}
	for(i = 0; i < 4; ++i) {
		memset(cs->cnt[i], 0, b->n * sizeof(int));
	}
	for(i = b->l; i < b->r; ++i) {
		char *seq = w->seqs[i].seq;
		for(j = 0; j < w->seqs[i].l_seq; ++j) {
			int pos = w->dis[i] + j;
			assert(pos < b->n);
			if(seq[j] > 3) continue;
			++cs->cnt[(int)seq[j]][pos];
		}
	}
	for(i = 0; i < 4; i++) {
		b->cnt[i] = malloc(b->n * sizeof(int));
	}
	for(i = 0; i < b->n; ++i) {
		int max = -1, max_id = 0;
		for(j = 0; j < 4; ++j) {
			b->cnt[j][i] = cs->cnt[j][i];
			if(cs->cnt[j][i] > max) {
				max = cs->cnt[j][i];
				max_id = j;
			}
		}
		b->s[i] = (uint8_t)max_id;
	}

	/* LAST-like */
	a->mem.n = 0;
	if(opt->max_mem_intv > 0) {
		int x = 0;
		while(x < b->n) {
			bwtintv_t m;
			x = bwt_seed_strategy1(bwt, b->n, b->s, x, opt->min_seed_len, opt->max_mem_intv, &m);
			if (m.x[2] > 0) {
				kv_push(bwtintv_t, a->mem, m);
			}
		}
	}

	/* SAL */
	mem_seed_v *seeds = &b->shared_seeds; kv_init(*seeds);
	for(i = 0; i < a->mem.n; ++i) {
		bwtintv_t *p = &a->mem.a[i];
		// append seeds of p to $seeds
		em2seeds(w->bwt, p, opt->max_occ, seeds);
	}
}

static void check_seeds_accuracy(const worker_t *w, int seq_id, int t_id, const mem_seed_v *seeds) {
	// Extract the simulated position from query name in format 'rname_read1Pos_read2Pos'
	const char *qname = w->seqs[seq_id].name;
	const int qlen = w->seqs[seq_id].l_seq;
	const bntseq_t *bns = w->bns;
	long i, j;
	char rname[64]; memset(rname, 0, sizeof(rname));
	for(j=0; qname[j] && qname[j]!='_'; j++) {
		rname[j] = qname[j];
	}
	if(!qname[j]) {
		fprintf(stderr, "Reference-name absent in wgsim read name.\n");
		abort();
	}
	int sim_pos1 = 0, sim_pos2 = 0;
	for(j=j+1; qname[j] && qname[j]!='_'; j++) {
		sim_pos1 *= 10; sim_pos1 += qname[j]-'0';
	}
	if(!qname[j]) {
		fprintf(stderr, "read1Pos absent in wgsim read name.\n");
		abort();
	}
	for(j=j+1; qname[j] && qname[j]!='_'; j++) {
		sim_pos2 *= 10; sim_pos2 += qname[j]-'0';
	}
	if(!qname[j]) {
		fprintf(stderr, "read2Pos absent in wgsim read name.\n");
		abort();
	}

	for (i = 0; i < seeds->n; i++) {
		mem_seed_t *p = &seeds->a[i];
		int rid = bns_intv2rid(bns, p->rbeg, p->rbeg + p->len);
		int rev_cs = w->is_rev[seq_id]; // Reversed in CS-generating.
		if(rid < 0) { continue; }
		bwtint_t pos;
		int is_rev;
		pos = bns_depos(bns, p->rbeg, &is_rev);
		if (is_rev) pos -= p->len - 1;
		is_rev = (rev_cs == 1) ?(is_rev ^ 1) :is_rev;
		pos = pos-bns->anns[rid].offset+1;
		// The seed is correct when close to simulated position enough.
		if(strcmp(bns->anns[rid].name, rname) != 0) continue;
		if(!is_rev && pos>=sim_pos1 && pos+p->len<sim_pos1+qlen) {
			zmp.correct_seeds_mt[t_id]++;
		}
		if( is_rev && pos>=sim_pos2 && pos+p->len<sim_pos2+qlen) {
			zmp.correct_seeds_mt[t_id]++;
		}
	}
}

void mez_collect_seeds(worker_t *w, int seq_id, int t_id) {
	const bwt_t *bwt = w->bwt;
	const mem_opt_t *opt = w->opt;
	uint8_t *seq = (uint8_t*)w->seqs[seq_id].seq;
	int_v *misbuf = &w->misbuf[t_id];
	int len = w->seqs[seq_id].l_seq;
	int cs_id = w->cs_id[seq_id];
	int dis = w->dis[seq_id];
	smem_aux_t *a = w->aux[t_id];

	/* 1. Seeding */
	int x, i, j, start_width = 1;
	int can_reuse = 1; // can reuse then zip-seeding, otherwise throw it into BWA-MEM seeding.
	misbuf->n = 0;
	a->mem.n = 0;

	if(cs_id == -1) {
		can_reuse = 0; // singletons
	} else {
		const r_batch_t *b = &w->batches.a[cs_id];
		if(b->shared_seeds.n != 0) {
			int last_err = -1, too_short = 1;
			for(i = 0; i < len; ++i) {
				if(b->s[dis + i] != seq[i]) {
					kv_push(int, *misbuf, i);
					if((i-1) - (last_err+1) + 1 >= opt->min_seed_len) {
						too_short = 0;
					}
					last_err = i;
				}
			}
			if((len-1) - (last_err+1) + 1 >= opt->min_seed_len) {
				too_short = 0;
			}
			int max_errors = (int)(len * 0.05 + 0.499);
			if(misbuf->n > max_errors || too_short == 1) { // mismatches too many or truncated segments too short.
				can_reuse = 0;
			}
		} else {
			can_reuse = 0; // don't have shared-seeds
		}
		/* Can't reuse shared-seeds, kick off it into BWA-MEM (this time cost should be a negligible part) */
	}

	/* Seeding */
	if(!can_reuse) {
		x = 0;
		while (x < len) {
			if (seq[x] < 4) {
				x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
				for (i = 0; i < a->mem1.n; ++i) {
					bwtintv_t *p = &a->mem1.a[i];
					if (intv_len(p) >= opt->min_seed_len) {
						kv_push(bwtintv_t, a->mem, *p);
					}
				}
			} else ++x;
		}
	} else {
		/* Reads in CS and not be kicked off, seeding on each mismatched position */
		const r_batch_t *b = &w->batches.a[cs_id];
		int last_x = -1;
		for(i = 0; i < misbuf->n; ++i) {
			int mis = misbuf->a[i];
			if(mis < last_x) continue; // This mismatch has been covered by previous EMs.
			if(seq[mis] > 3) continue; // Can't seeding on an ambiguous base N
			int sum = b->cnt[0][dis+mis] + b->cnt[1][dis+mis] + b->cnt[2][dis+mis] + b->cnt[3][dis+mis];
			if(b->cnt[seq[mis]][dis+mis] <= sum / 4) continue; // The mismatched base takes few percentage.
			last_x = bwt_smem1(bwt, len, seq, mis, start_width, &a->mem1, a->tmpv);
			for(j = 0; j < a->mem1.n; ++j) {
				bwtintv_t *p = &a->mem1.a[j];
				if(intv_len(p) >= opt->min_seed_len) {
					kv_push(bwtintv_t, a->mem, *p);
				}
			}
		}
	}

	/* 2. Re-seeding */
	int old_n = (int)a->mem.n;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	for(i = 0; i < old_n; ++i) {
		bwtintv_t *p = &a->mem.a[i];
		int start = (int)(p->info>>32), end = (int32_t)p->info;
		if (end - start < split_len || p->x[2] > opt->split_width) continue;
		bwt_smem1(bwt, len, seq, (start + end)>>1, p->x[2] + 1, &a->mem1, a->tmpv);
		for(j = 0; j < a->mem1.n; ++j) {
			if(intv_len(&a->mem1.a[j]) >= opt->min_seed_len) {
				kv_push(bwtintv_t, a->mem, a->mem1.a[j]);
			}
		}
	}

	/* 3. LAST-like */
	if(can_reuse == 0 && opt->max_mem_intv > 0) {
		x = 0;
		while(x < len) {
			if(seq[x] < 4) {
				bwtintv_t m;
				x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
				if(m.x[2] > 0) {
					kv_push(bwtintv_t, a->mem, m);
				}
			} else {
				++x;
			}
		}
	}

	if(can_reuse == 0) {
		w->cs_id[seq_id] = -1; // Mark this for seeds sorting.
		ks_introsort(mez_intv, a->mem.n, a->mem.a);
	}

	/* 4. SAL */
	mem_seed_v *seeds = &w->seq_seeds[t_id]; seeds->n = 0;
	for (i = 0; i < a->mem.n; ++i) {
		bwtintv_t *p = &a->mem.a[i];
		em2seeds(bwt, p, opt->max_occ, seeds);
	}

	if(can_reuse == 0) {
		check_seeds_accuracy(w, seq_id, t_id, seeds);
		zmp.seeds_mt[t_id] += seeds->n;
		return ;
	}

	/* 5. Reuse shared seeds in CS */
	const r_batch_t *b = &w->batches.a[cs_id];
	const mem_seed_v *ss = &b->shared_seeds;
	for(i = 0; i < ss->n; ++i) { // seeds with query-begin < dis + len
		mem_seed_t *s = &ss->a[i];
		if(s->qbeg >= dis + len || s->qbeg + s->len <= dis) { // has no intersection
			continue;
		}
		int start = s->qbeg - dis, end = s->qbeg + s->len - dis; // [start, end) had fall down query read
		int cut_l = (start > 0) ?start :0; // cut [start, end)
		int cut_r = (end < len) ?end :len;
		assert(cut_l < cut_r);
		if(cut_r - cut_l < opt->min_seed_len) continue;
		int_v *misv = &w->seedmis[t_id]; misv->n = 0;
		kv_push(int, *misv, cut_l - 1);
		for(j = 0; j < misbuf->n; ++j) {
			int mis = misbuf->a[j];
			if(mis >= cut_l && mis < cut_r) {
				kv_push(int, *misv, mis);
			}
		}
		kv_push(int, *misv, cut_r);
		for(j = 0; j < misv->n - 1; ++j) {
			int qb = misv->a[j] + 1, qe = misv->a[j+1]; // [qb, qe)
			if(qe - qb < opt->min_seed_len) continue;
			mem_seed_t rs;
			rs.rbeg = s->rbeg + (qb - start);
			rs.qbeg = qb;
			rs.len = qe - qb;
			rs.score = rs.len; // reused seeds didn't call em2seeds, so we need to assign score manually.
			kv_push(mem_seed_t, *seeds, rs);
		}
	}
	check_seeds_accuracy(w, seq_id, t_id, seeds);
	zmp.seeds_mt[t_id] += seeds->n;
}

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(mezchn, mem_chain_t, chain_cmp)

mem_chain_v mez_chain(worker_t *w, int seq_id, int t_id) {
	mem_seed_v *seeds = &w->seq_seeds[t_id];
	int len = w->seqs[seq_id].l_seq;
	const bntseq_t *bns = w->bns;
	const mem_opt_t *opt = w->opt;
	const smem_aux_t *aux = w->aux[t_id];

	mem_chain_v chn; kv_init(chn);
	if(seeds->n == 0) {
		return chn;
	}
	if(w->cs_id[seq_id] != -1) { // Singletons and kicked-off reads are sorted by EM:{qb, qe}
		ks_introsort(mem_seeds, seeds->n, seeds->a);
	}
	kbtree_t(mezchn) *tree = kb_init(mezchn, KB_DEFAULT_SIZE);
	int i, b, e, l_rep = 0;

	/* compute frac_rep */
	for (i = 0, b = e = l_rep = 0; i < aux->mem.n; ++i) {
		bwtintv_t *p = &aux->mem.a[i];
		int sb = (p->info>>32), se = (uint32_t)p->info;
		if (p->x[2] <= opt->max_occ) continue;
		if (sb > e) l_rep += e - b, b = sb, e = se;
		else e = e > se? e : se;
	}
	l_rep += e - b;

	for(i = 0; i < seeds->n; ++i) {
		const mem_seed_t *s = &seeds->a[i];
		mem_chain_t tmp, *lower, *upper;
		tmp.pos = s->rbeg;
		int rid, to_add = 0;
		rid = bns_intv2rid(bns, s->rbeg, s->rbeg + s->len);
		if(rid < 0) continue;
		if(kb_size(tree)) {
			kb_intervalp(mezchn, tree, &tmp, &lower, &upper);
			if(!lower || !test_and_merge(opt, bns->l_pac, lower, s, rid)) {
				to_add = 1;
			}
		} else {
			to_add = 1;
		}
		if(to_add) {
			tmp.n = 1; tmp.m = 4;
			tmp.seeds = calloc(tmp.m, sizeof(mem_seed_t));
			tmp.seeds[0] = *s;
			tmp.rid = rid;
			tmp.is_alt = !!bns->anns[rid].is_alt;
			kb_putp(mezchn, tree, &tmp);
		}
	}

	kv_resize(mem_chain_t, chn, kb_size(tree));

#define traverse_func(p_) (chn.a[chn.n++] = *(p_))
	__kb_traverse(mem_chain_t, tree, traverse_func);
#undef traverse_func

	for(i = 0; i < chn.n; ++i) {
		chn.a[i].frac_rep = (float)l_rep / len;
	}
	if(bwa_verbose >= 4) {
		printf("* fraction of repetitive seeds: %.3f\n", (float)l_rep / len);
	}
	kb_destroy(chn, tree);

	return chn;
}

mem_alnreg_v mez_align1_core(worker_t *w, int seq_id, int t_id) {
	const mem_opt_t *opt = w->opt;
	const bntseq_t *bns = w->bns;
	const uint8_t *pac = w->pac;
	uint8_t *seq = (uint8_t*)w->seqs[seq_id].seq;
	int l_seq = w->seqs[seq_id].l_seq;

	/* Collecting seeds. */
	mez_collect_seeds(w, seq_id, t_id);

	/* chaining collected seeds */
	mem_chain_v chn = mez_chain(w, seq_id, t_id);

	/* chains filtering */
	chn.n = (size_t)mem_chain_flt(opt, chn.n, chn.a);

	/* Chains to aligned regions */
	zmp.chains_mt[t_id] += chn.n;
	mem_alnreg_v regs; kv_init(regs);
	int i;
	for (i = 0; i < chn.n; ++i) {
		mem_chain_t *p = &chn.a[i];
		if (bwa_verbose >= 4) err_printf("* ---> Processing chain(%d) <---\n", i);
		mem_chain2aln(opt, bns, pac, l_seq, seq, p, &regs);
		free(p->seeds);
	}
	free(chn.a);
	zmp.regs_mt[t_id] += regs.n;

	/* regions patch */
	regs.n = (size_t)mem_sort_dedup_patch(opt, bns, pac, seq, regs.n, regs.a);
	zmp.patch_regs_mt[t_id] += regs.n;

	if (bwa_verbose >= 4) {
		err_printf("* %ld chains remain after removing duplicated chains\n", regs.n);
		for (i = 0; i < regs.n; ++i) {
			mem_alnreg_t *p = &regs.a[i];
			printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
		}
	}
	for (i = 0; i < regs.n; ++i) {
		mem_alnreg_t *p = &regs.a[i];
		if (p->rid >= 0 && bns->anns[p->rid].is_alt)
			p->is_alt = 1;
	}
	return regs;
}

static inline int get_rlen(int n_cigar, const uint32_t *cigar)
{
	int k, l;
	for (k = l = 0; k < n_cigar; ++k) {
		int op = cigar[k]&0xf;
		if (op == 0 || op == 2)
			l += cigar[k]>>4;
	}
	return l;
}

static inline void add_cigar(const mem_opt_t *opt, mem_aln_t *p, kstring_t *str, int which)
{
	int i;
	if (p->n_cigar) { // aligned
		for (i = 0; i < p->n_cigar; ++i) {
			int c = p->cigar[i]&0xf;
			if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
				c = which? 4 : 3; // use hard clipping for supplementary alignments
			kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
		}
	} else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
}

void mez_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_, uint8_t is_rev)
{
	int i, l_name;
	mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

	if (m_) mtmp = *m_, m = &mtmp;
	// set flag
	p->flag |= m? 0x1 : 0; // is paired in sequencing
	p->flag |= p->rid < 0? 0x4 : 0; // is mapped
	p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
	if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
		p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
	if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
		m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
	p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
	p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand
	// print up to CIGAR
	l_name = strlen(s->name);
	ks_resize(str, str->l + s->l_seq + l_name + (s->qual? s->l_seq : 0) + 20);
	kputsn(s->name, l_name, str); kputc('\t', str); // QNAME
	// We might have reversed the read in CS-generating stage.
	p->flag ^= is_rev ?0x10 :0;
	kputw((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc('\t', str); // FLAG
	// Undo the modification on p->flag to avoid the impact on SEQ and QAUL, only changing the RC flag in SAM is ok and safe.
	p->flag ^= is_rev ?0x10 :0;
	if (p->rid >= 0) { // with coordinate
		kputs(bns->anns[p->rid].name, str); kputc('\t', str); // RNAME
		kputl(p->pos + 1, str); kputc('\t', str); // POS
		kputw(p->mapq, str); kputc('\t', str); // MAPQ
		add_cigar(opt, p, str, which);
	} else kputsn("*\t0\t0\t*", 7, str); // without coordinte
	kputc('\t', str);

	// print the mate position if applicable
	if (m && m->rid >= 0) {
		if (p->rid == m->rid) kputc('=', str);
		else kputs(bns->anns[m->rid].name, str);
		kputc('\t', str);
		kputl(m->pos + 1, str); kputc('\t', str);
		/* infer insert size */
		if (p->rid == m->rid) {
			int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
			int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
			if (m->n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
			else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
		} else kputc('0', str);
	} else kputsn("*\t0\t0", 5, str);
	kputc('\t', str);

	// print SEQ and QUAL
	if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
		kputsn("*\t*", 3, str);
	} else if (!p->is_rev) { // the forward strand
		int qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // have cigar && not the primary alignment && not softclip all
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	} else { // the reverse strand
		int qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	}

	// print optional tags
	if (p->n_cigar) {
		kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
		kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
	}
	if (m && m->n_cigar) { kputsn("\tMC:Z:", 6, str); add_cigar(opt, m, str, which); }
	if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
	if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(p->sub, str); }
	if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
	if (!(p->flag & 0x100)) { // not multi-hit
		for (i = 0; i < n; ++i)
			if (i != which && !(list[i].flag&0x100)) break;
		if (i < n) { // there are other primary hits; output them
			kputsn("\tSA:Z:", 6, str); // SA: Supplementary Alignment.
			for (i = 0; i < n; ++i) {
				const mem_aln_t *r = &list[i];
				int k;
				if (i == which || (r->flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hi
				kputs(bns->anns[r->rid].name, str); kputc(',', str);
				kputl(r->pos+1, str); kputc(',', str);
				kputc("+-"[r->is_rev], str); kputc(',', str);
				for (k = 0; k < r->n_cigar; ++k) {
					kputw(r->cigar[k]>>4, str); kputc("MIDSH"[r->cigar[k]&0xf], str);
				}
				kputc(',', str); kputw(r->mapq, str);
				kputc(',', str); kputw(r->NM, str);
				kputc(';', str);
			}
		}
		if (p->alt_sc > 0)
			ksprintf(str, "\tpa:f:%.3f", (double)p->score / p->alt_sc);
	}
	if (p->XA) { kputsn("\tXA:Z:", 6, str); kputs(p->XA, str); }
	if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
	if ((opt->flag&MEM_F_REF_HDR) && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
		int tmp;
		kputsn("\tXR:Z:", 6, str);
		tmp = str->l;
		kputs(bns->anns[p->rid].anno, str);
		for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
			if (str->s[i] == '\t') str->s[i] = ' ';
	}
	kputc('\n', str);
}

void mez_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s,
                 mem_alnreg_v *a, int extra_flag, const mem_aln_t *m, uint8_t is_rev, int t_id)
{
	kstring_t str;
	kvec_t(mem_aln_t) aa;
	int k, l;
	char **XA = 0;

	if (!(opt->flag & MEM_F_ALL)) {
		XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
	}
	kv_init(aa);
	str.l = str.m = 0; str.s = 0;

	for (k = l = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		mem_aln_t *q;
		if (p->score < opt->T) continue;

		if (p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))) continue;

		if (p->secondary >= 0 && p->secondary < INT_MAX && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;

		q = kv_pushp(mem_aln_t, aa);
		*q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
		zmp.reg2aln_mt[t_id]++;

		assert(q->rid >= 0); // this should not happen with the new code

		q->XA = XA? XA[k] : 0;

		q->flag |= extra_flag; // flag secondary

		if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score

		/* Supplementary alignments. primary hits in chimeric alignment. */
		if (l && p->secondary < 0) // if supplementary
			q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;

		if (!(opt->flag & MEM_F_KEEP_SUPP_MAPQ) && l && !p->is_alt && q->mapq > aa.a[0].mapq)
			q->mapq = aa.a[0].mapq; // lower mapq for supplementary mappings, unless -5 or -q is applied
		++l;
	}
	if (aa.n == 0) { // no alignments good enough; then write an unaligned record
		mem_aln_t t;
		t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
		t.flag |= extra_flag;
		mem_aln2sam(opt, bns, &str, s, 1, &t, 0, m);
	} else {
		for (k = 0; k < aa.n; ++k)
			mez_aln2sam(opt, bns, &str, s, aa.n, aa.a, k, m, is_rev);
		for (k = 0; k < aa.n; ++k) free(aa.a[k].cigar);
		free(aa.a);
	}
	s->sam = str.s;
	if (XA) {
		for (k = 0; k < a->n; ++k) free(XA[k]);
		free(XA);
	}
}

static void mez_worker1(void *data, int i, int tid) {
	worker_t *w = (worker_t*)data;
	if (!(w->opt->flag&MEM_F_PE)) {
		if (bwa_verbose >= 4) printf("=====> Processing read '%s' <=====\n", w->seqs[i].name);
		mem_alnreg_v regs = mez_align1_core(w, i, tid);

		/* SAM-format */
		if (bwa_verbose >= 4) printf("=====> Finalizing read '%s' <=====\n", w->seqs[i].name);
		mem_mark_primary_se(w->opt, regs.n, regs.a, w->n_processed + i);
		if (w->opt->flag & MEM_F_PRIMARY5) mem_reorder_primary5(w->opt->T, &regs);
		mez_reg2sam(w->opt, w->bns, w->pac, &w->seqs[i], &regs, 0, 0, w->is_rev[i], tid);
		free(regs.a);
	}
}


void zip_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns,
                      const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs)
{
	extern void kt_for(int n_threads, void (*func)(void*,int,int), void *data, int n);
	worker_t w;
	mem_pestat_t pes[4];
	double ctime, rtime;
	int i;

	ctime = cputime(); rtime = realtime();
	w.opt = opt; w.bwt = bwt; w.bns = bns; w.pac = pac;
	w.seqs = seqs; w.n_processed = n_processed;
	w.n = n;

	/* Get consensus sequences, parallel process each chunk */
	w.is_rev = calloc(n, sizeof(uint8_t));
	w.dis = malloc(n * sizeof(int)); memset(w.dis, -1, n * sizeof(int ));

	int_v input_boundary; kv_init(input_boundary); kv_push(int, input_boundary, 0);
	size_t bytes = 0;
	int cnt_reads = 0;
	for(i = 0; i < n; ++i) {
		++cnt_reads;
		bytes += seqs[i].l_seq;
		if(i == n - 1 || (bytes >= opt->chunk_size && (cnt_reads & 1) == 0)) {
			kv_push(int, input_boundary, i + 1);
			cnt_reads = 0;
			bytes = 0;
		}
	}
	assert(input_boundary.n-1 <= w.opt->n_threads);
	kv_init(w.chunk_intv); kv_push(int, w.chunk_intv, 0);
	for (i = 0; i < input_boundary.n - 1; i++) {
		int bdy_l = input_boundary.a[i];
		int bdy_r = input_boundary.a[i+1];
		int j;
		bytes = 0;
		for (j = bdy_l; j < bdy_r; j++) {
			bytes += seqs[j].l_seq;
			if(j == bdy_r-1 || bytes >= opt->chunk_size / 10) {
				kv_push(int, w.chunk_intv, j + 1);
				bytes = 0;
			}
		}
	}
	free(input_boundary.a);

	w.back_batches = calloc(opt->n_threads, sizeof(r_batch_v));
	kt_for(opt->n_threads, gen_cs, &w, w.chunk_intv.n-1);
	free(w.chunk_intv.a);

	// Sum up CSs in all threads.
	kv_init(w.batches);
	for(i = 0; i < opt->n_threads; ++i) {
		const r_batch_v *p = &w.back_batches[i];
		int j;
		for(j = 0; j < p->n; ++j) {
			kv_push(r_batch_t, w.batches, p->a[j]);
		}
		free(p->a);
	}
	free(w.back_batches);
	double gencs_c = cputime(), gencs_r = realtime();
	zmp.t_gencs[0] += gencs_c - ctime; zmp.t_gencs[1] += gencs_r - rtime;

	/* Seeding on Consensus Sequences */
	w.aux = malloc(opt->n_threads * sizeof(smem_aux_t*));  for (i = 0; i < opt->n_threads; i++) w.aux[i] = smem_aux_init();
	w.cs_aux = malloc(opt->n_threads * sizeof(cs_aux_t*)); for (i = 0; i < opt->n_threads; i++) w.cs_aux[i] = cs_aux_init();
	w.cs_id = malloc(n * sizeof(int)); memset(w.cs_id, -1, n * sizeof(int));
	kt_for(opt->n_threads, cs_seeding, &w, (int) w.batches.n);
	double csseed_c = cputime(), csseed_r = realtime();
	zmp.t_csseed[0] += csseed_c - gencs_c; zmp.t_csseed[1] += csseed_r - gencs_r;

	/* Collect profile */
	do {
		zmp.reads_n += n;
		zmp.read_len = seqs[0].l_seq;
		zmp.singletons += n;
		for (i = 0; i < w.batches.n; i++) {
			zmp.singletons -= (w.batches.a[i].r - w.batches.a[i].l);
			zmp.cs_bases += w.batches.a[i].n;
		}
	} while (0);

	/* Reusing and supplementary seeding (process singletons) */
	w.misbuf = calloc(opt->n_threads, sizeof(int_v));
	w.seedmis = calloc(opt->n_threads, sizeof(int_v));
	w.seq_seeds = calloc(opt->n_threads, sizeof(mem_seed_v));

	kt_for(opt->n_threads, mez_worker1, &w, (opt->flag & MEM_F_PE) ? n >> 1 : n); // find mapping positions

	// Temp memory free
	free(w.dis);
	for (i = 0; i < w.batches.n; i++)  {
		const r_batch_t *p = &w.batches.a[i];
		free(p->s);
		free(p->shared_seeds.a);
		free(p->cnt[0]); free(p->cnt[1]); free(p->cnt[2]); free(p->cnt[3]);
	}
	free(w.batches.a);
	for (i = 0; i < opt->n_threads; i++) { smem_aux_destroy(w.aux[i]); } free(w.aux);
	for (i = 0; i < opt->n_threads; i++) { cs_aux_destroy(w.cs_aux[i]); } free(w.cs_aux);
	free(w.cs_id);
	for (i = 0; i < opt->n_threads; i++) { free(w.misbuf[i].a); } free(w.misbuf);
	for (i = 0; i < opt->n_threads; i++) { free(w.seedmis[i].a); } free(w.seedmis);
	for (i = 0; i < opt->n_threads; i++) { free(w.seq_seeds[i].a); }  free(w.seq_seeds);

	double worker1_c = cputime(), worker1_r = realtime();
	zmp.t_worker1[0] += worker1_c - csseed_c; zmp.t_worker1[1] += worker1_r - csseed_r;

	free(w.is_rev);

	fprintf(stderr, "[M::%s] Processed %d reads in %.2f CPU seconds, %.2f seconds elapsed.\n", __func__ , n, cputime()-ctime, realtime()-rtime);
}
