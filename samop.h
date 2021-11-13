//
// Created by ixiaohu on 2021/11/13.
//

#ifndef ZIP_SEEDING_SAMOP_H
#define ZIP_SEEDING_SAMOP_H

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

#ifdef __cplusplus
extern "C" {
#endif
	sam_line_t fetch_samline1(gzFile f);
#ifdef __cplusplus
}
#endif

#endif //ZIP_SEEDING_SAMOP_H
