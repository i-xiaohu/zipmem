//
// Created by ixiaohu on 2021/10/14.
//

#include <stdio.h>

#include "cstl/kvec.h"
#include "bwalib/bwa.h"
#include "bwalib/kseq.h"

KSEQ_DECLARE(gzFile)

void view_seeds(const char *idx_name, const char *fastq_fn, const char *seed_fn) {
	bwaidx_t *idx = bwa_idx_load_from_disk(idx_name, BWA_IDX_ALL);
	gzFile fp = gzopen(fastq_fn, "r");
	kseq_t *ks = kseq_init(fp);
	FILE *fseeds = fopen(seed_fn, "r");
	long reads_n = 0;
	long seeds_n = 0, no_mismatch = 0, mismatches_n = 0;
	while (1) {
		int i, j, k, h, n_seqs;
		long size = 0;
		bseq1_t *seqs = bseq_read(240 * 1000 * 1000, &n_seqs, ks, NULL);
		if (seqs == NULL) break;
		reads_n += n_seqs;
		for (i = 0; i < n_seqs; ++i) size += seqs[i].l_seq;
		long total_bytes = 0;
		for (i = 0; i < n_seqs; i++) {
			bseq1_t *read = &seqs[i];
			for (j = 0; j < read->l_seq; j++) read->seq[j] = nst_nt4_table[read->seq[j]];
			long seed_bytes = 0;
			fread(&seed_bytes, sizeof(long), 1, fseeds);
			total_bytes += seed_bytes;
			uint8_t *seeds = malloc(seed_bytes);
			fread(seeds, sizeof(uint8_t), seed_bytes, fseeds);
			long pointer = 0;
			int mem_n = *(int*)(seeds + pointer); pointer += sizeof(int);
			for (j = 0; j < mem_n; j++) {
				int sb = *(int*)(seeds + pointer); pointer += sizeof(int);
				int sl = *(int*)(seeds + pointer); pointer += sizeof(int);
				long sa_size = *(long*)(seeds + pointer); pointer += sizeof(long);
				int true_occ = sa_size < 500 ?sa_size :500;
				seeds_n += true_occ;
				for (k = 0; k < true_occ; k++) {
					long rb = *(long*)(seeds + pointer); pointer += sizeof(long);
					int rid = bns_intv2rid(idx->bns, rb, rb + sl);
					if (rid < 0) continue;
					int64_t get_len = 0;
					uint8_t *ref = bns_get_seq(idx->bns->l_pac, idx->pac, rb, rb + sl, &get_len);
					if (get_len != sl) {
						fprintf(stderr, "%d %ld\n", sl, get_len);
					}
					assert(get_len == sl);
					int mis = 0;
					for (h = 0; h < sl; h++) {
						if (read->seq[sb + h] != ref[h]) mis++;
					}
					free(ref);
					if (!mis) no_mismatch++;
					else mismatches_n += mis;
				}
			}
			if (pointer != seed_bytes) {
				fprintf(stderr, "%ld %ld\n", pointer, seed_bytes);
			}
			assert(pointer == seed_bytes);
			free(seeds);
			free(read->name); free(read->comment); free(read->seq); free(read->qual);
		}
		free(seqs);
		fprintf(stderr, "[%s] Processed %d reads (%ld Mbp) and %ld MB seeds\n",
		        __func__, n_seqs, (long)size/1000/1000, total_bytes/1024/1024);
	}
	bwa_idx_destroy(idx);
	kseq_destroy(ks);
	gzclose(fp);
	fclose(fseeds);
	fprintf(stderr, "reads number\t%ld\n", reads_n);
	fprintf(stderr, "seeds number\t%.2f\n", 1.0 * seeds_n / reads_n);
	fprintf(stderr, "Exactly matched seeds\t%.2f\n", 100.0 * no_mismatch / seeds_n);
	fprintf(stderr, "Mismatched number of noisy seeds\t%.2f\n", 1.0 * mismatches_n / (seeds_n-no_mismatch));
}

static int usage() {
	fprintf(stderr, "\n");
	fprintf(stderr, "Check seeds and SAM alignments of zipmem/bwamem\n");
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "    check seed <zipmem.seeds> <bwamem.seeds>\n");
	fprintf(stderr, "    check sam <zipmem.sam> <bwamem.sam>\n");
	fprintf(stderr, "\n");
	return 0;
}

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

void check_seeds(const char *zip_fn, const char *mem_fn) {

}

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

static void output_sam(const sam_line_t *s) {
	int i;
	fprintf(stderr, "QName: %s\n", s->qname);
	fprintf(stderr, "Flag:  ");
	for (i = 11; i >= 0; i--)
		if (s->flag & (1<<i)) fprintf(stderr, "1");
		else fprintf(stderr, "0");
	fprintf(stderr, "\n");
	fprintf(stderr, "RName: %s\n", s->rname);
	fprintf(stderr, "Pos:   %d\n", s->pos);
	fprintf(stderr, "MapQ:  %d\n", s->mapq);
	fprintf(stderr, "CIGAR: %s\n", s->cigar);
	fprintf(stderr, "SEQ:   %s\n", s->seq);
	fprintf(stderr, "NM:    %d\n", s->nm);
	fprintf(stderr, "AS:    %d\n", s->as);
	fprintf(stderr, "\n");
}

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

void check_sam(const char *zip_fn, const char *mem_fn) {
	gzFile fzip = gzopen(zip_fn, "r"); assert(fzip != NULL);
	gzFile fmem = gzopen(mem_fn, "r"); assert(fmem != NULL);
	// Only consider primary alignments
	long pos_eq = 0, primary_n = 0;
	long as_eq = 0;
	long nm_eq = 0;
	while (1) {
		sam_line_t zip_sam = fetch_samline1(fzip);
		if (gzeof(fzip)) break;
		if (zip_sam.data == NULL) continue;
		if (sec_ali(zip_sam.flag) || sup_ali(zip_sam.flag)) {
			free(zip_sam.data);
			continue;
		}
		sam_line_t mem_sam;
		while (1) {
			mem_sam = fetch_samline1(fmem);
			if (gzeof(fmem)) break;
			if (mem_sam.data == NULL) continue;
			if (sec_ali(mem_sam.flag) || sup_ali(mem_sam.flag)) {
				free(mem_sam.data);
				continue;
			}
			break;
		}
		if (gzeof(fmem)) break;
		assert(strcmp(zip_sam.qname, mem_sam.qname) == 0);
		primary_n++;
		if (!strcmp(zip_sam.rname, mem_sam.rname) && zip_sam.pos == mem_sam.pos) {
			pos_eq++;
		}
		if (zip_sam.as == mem_sam.as) as_eq++;
		if (zip_sam.nm == mem_sam.nm) nm_eq++;
		free(zip_sam.data); free(mem_sam.data);
	}
	gzclose(fzip); gzclose(fmem);
	fprintf(stderr, "Primary: %ld\n", primary_n);
	fprintf(stderr, "POS:\t%.2f %% = %ld / #pri\n", 100.0 * pos_eq / primary_n, pos_eq);
	fprintf(stderr, "AS: \t%.2f %% = %ld / #pri\n", 100.0 * as_eq / primary_n, as_eq);
	fprintf(stderr, "NM: \t%.2f %% = %ld / #pri\n", 100.0 * nm_eq / primary_n, nm_eq);
}

typedef kvec_t(bseq1_t) bseq1_v;
void view_sam(const char *ref_fn, const char *zip_fn) {
	gzFile fzip = gzopen(zip_fn, "r"); assert(fzip != NULL);
	gzFile fp = gzopen(ref_fn, "r");
	kseq_t *ks = kseq_init(fp);

	bseq1_v seqs; kv_init(seqs);
	while (kseq_read(ks) >= 0) {
		bseq1_t seq;
		seq.name = strdup(ks->name.s);
		seq.seq = strdup(ks->seq.s);
		kv_push(bseq1_t, seqs, seq);
	}
	kseq_destroy(ks); gzclose(fp);
	fprintf(stderr, "Input %ld reference sequences\n", seqs.n);

	int i, j;
	while (1) {
		sam_line_t zip_sam = fetch_samline1(fzip);
		if (gzeof(fzip)) break;
		if (zip_sam.data == NULL) continue;
		if (sec_ali(zip_sam.flag) || sup_ali(zip_sam.flag)) {
			free(zip_sam.data);
			continue;
		}

		char *que = zip_sam.seq;
		char *ref = NULL;
		for (i = 0; i < seqs.n; i++) {
			bseq1_t *seq = &seqs.a[i];
			if (!strcmp(seq->name, zip_sam.rname)) {
				ref = seq->seq;
				break;
			}
		}
		assert(ref != NULL);

		int match = 1, mis_pen = 4;
		int o_del = 6, e_del = 1;
		int o_ins = 6, e_ins = 1;
		int clip_pen = 5;
		int temp = 0, pr = zip_sam.pos-1, pq = 0, score = 0;
		for (i = 0; zip_sam.cigar[i]; i++) {
			char c = zip_sam.cigar[i];
			if (c >= '0' && c <= '9') {
				temp *= 10;
				temp += c - '0';
			} else {
				if (c == 'M') {
					for (j = 0; j < temp; j++) {
						if (que[pq] == 'N') score -= 1;
						else if (que[pq] != ref[pr]) score -= mis_pen;
						else score += match;
						pq++; pr++;
					}
				} else if (c == 'I') {
					score -= o_ins + temp * e_ins;
					pq += temp;
				} else if (c == 'D') {
					score -= o_del + temp * e_del;
					pr += temp;
				} else if (c == 'S') {
					score -= clip_pen;
					pq += temp;
				} else abort();
				temp = 0;
			}
		}
		assert(pq == strlen(zip_sam.seq));
		if (score != zip_sam.as) { // This is actually not identical to BWA-MEM
			output_sam(&zip_sam);
			fprintf(stderr, "Score I calculated: %d\n", score);
			for (i = zip_sam.pos-1; i < pr; i++) {
				fprintf(stderr, "%c", ref[i]);
			}
			fprintf(stderr, "\n");
			for (i = zip_sam.pos-1; i < pr; i++) {
				if (ref[i] != que[i-zip_sam.pos+1]) fprintf(stderr, "-");
				else fprintf(stderr, "*");
			}
			fprintf(stderr, "\n");
			break;
		}
		free(zip_sam.data);
	}
	gzclose(fzip);
}

int main(int argc, char *argv[]) {
	if (argc == 1) return usage();
	if (!strcmp(argv[1], "seed") && argc == 4) check_seeds(argv[2], argv[3]);
	else if (!strcmp(argv[1], "seed") && argc == 5) view_seeds(argv[2], argv[3], argv[4]);
	else if (!strcmp(argv[1], "2sam")) check_sam(argv[2], argv[3]);
	else if (!strcmp(argv[1], "1sam")) view_sam(argv[2], argv[3]);
	else return usage();
	return 0;
}
