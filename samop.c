//
// Created by ixiaohu on 2021/11/13.
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#include "samop.h"

sam_line_t fetch_samline1(gzFile f) {
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