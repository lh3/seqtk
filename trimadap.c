#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <zlib.h>
#include "ksw.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct {
	int type, len;
	uint8_t *seq;
	kswq_t *qp;
	uint64_t cnt;
} ta_adap_t;

int main(int argc, char *argv[])
{
	int n_adaps, m_adaps;
	int c, i, j, k, from_stdin;
	int sa = 1, sb = 2, go = 1, ge = 3, type = 1;
	int min_sc = 15, min_len = 10;
	double max_diff = 0.15;
	ta_adap_t *adaps;
	kseq_t *ks;
	gzFile fp;
	int8_t mat[25];
	kstring_t str = {0,0,0};

	n_adaps = m_adaps = 0; adaps = 0;
	while ((c = getopt(argc, argv, "5:3:s:t:l:")) >= 0) {
		if (c == '5' || c == '3') {
			ta_adap_t *p;
			if (m_adaps == n_adaps) {
				m_adaps = m_adaps? m_adaps<<1 : 4;
				adaps = realloc(adaps, m_adaps * sizeof(ta_adap_t));
			}
			p = &adaps[n_adaps++];
			p->seq = (uint8_t*)strdup(optarg);
			p->type = c - '0';
		} else if (c == 't') {
			if (strcmp(optarg, "ilpe") == 0) type = 1;
		} else if (c == 's') min_sc = atoi(optarg);
		else if (c == 'd') max_diff = atof(optarg);
		else if (c == 'l') min_len = atoi(optarg);
	}

	// preset
	if (type == 1 && n_adaps == 0) {
		m_adaps = n_adaps = 3;
		adaps = malloc(m_adaps * sizeof(ta_adap_t));
		adaps[0].type = 5; adaps[0].seq = (uint8_t*)strdup("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT");
		adaps[1].type = 3; adaps[1].seq = (uint8_t*)strdup("AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC");
		adaps[2].type = 3; adaps[2].seq = (uint8_t*)strdup("ATCTCGTATGCCGTCTTCTGCTTG");
	}

	// update adapter info
	for (j = 0; j < n_adaps; ++j) {
		ta_adap_t *p = &adaps[j];
		p->len = strlen((char*)p->seq);
		p->qp = 0;
		p->cnt = 0;
		for (i = 0; i < p->len; ++i)
			p->seq[i] = seq_nt4_table[(uint8_t)p->seq[i]];
	}

	from_stdin = !isatty(fileno(stdin));
	if (optind == argc && !from_stdin) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   trimadap [options] <in.fq>\n\n");
		fprintf(stderr, "Options: -5 STR     5'-end adapter\n");
		fprintf(stderr, "         -3 STR     3'-end adapter\n");
		fprintf(stderr, "         -l INT     min length [%d]\n", min_len);
		fprintf(stderr, "         -s INT     min score [%d]\n", min_sc);
		fprintf(stderr, "         -d FLOAT   max difference [%.3f]\n", max_diff);
		fprintf(stderr, "\n");
		return 1; // FIXME: memory leak
	}

	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? sa : -sb;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;

	fp = optind < argc && strcmp(argv[optind], "-")? gzopen(argv[optind], "rb") : gzdopen(fileno(stdin), "rb");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		if (str.m < ks->seq.m) {
			str.m = ks->seq.m;
			str.s = realloc(str.s, str.m);
		}
		str.l = ks->seq.l;
		for (i = 0; i < ks->seq.l; ++i)
			str.s[i] = seq_nt4_table[(uint8_t)ks->seq.s[i]];
		for (j = 0; j < n_adaps; ++j) {
			kswr_t r;
			double diff;
			int type;
			ta_adap_t *p = &adaps[j];
			r = ksw_align(p->len, p->seq, str.l, (uint8_t*)str.s, 5, mat, go, ge, KSW_XBYTE|KSW_XSTART|(min_len*sa), &p->qp);
			++r.te; ++r.qe; // change to 0-based
			k = r.qe - r.qb < r.te - r.tb? r.qe - r.qb : r.te - r.tb;
			diff = (double)(k * sa - r.score) / sb / k;
			//printf("%d:%.3f [%d,%d):%d <=> [%d,%d):%d\n", r.score, diff, r.qb, r.qe, p->len, r.tb, r.te, (int)str.l);
			if (r.qb <= r.tb && p->len - r.qe <= str.l - r.te) { // contained
				if (r.qb * sa > sa + sb) continue;
				if ((p->len - r.qe) * sa > sa + sb) continue;
				type = 1;
			} else if (r.qb <= r.tb) { // 3' overlap
				if (r.qb * sa > sa + sb) continue;
				if ((str.l - r.te) * sa > sa + sb) continue;
				type = 2;
			} else {
				if ((p->len - r.qe) * sa > sa + sb) continue;
				if (r.tb * sa > sa + sb) continue;
				type = 3;
			}
			if (p->type == 5) {
				if (r.tb == 0 && r.qe == p->len && (r.te - r.tb) * sa == r.score)
					type = 4;
			} else if (p->type == 3) {
				if (r.qb == 0 && r.te == str.l && (r.te - r.tb) * sa == r.score)
					type = 4;
			}
			if (type == 4) {
				if (r.te - r.tb < min_len) continue;
			} else {
				if (r.score < min_sc || diff > max_diff) continue;
			}
			++p->cnt;
			if (p->type == 5) {
				k = r.te + (p->len - r.qe);
				k = k < str.l? k : str.l;
				for (i = 0; i < k; ++i) ks->seq.s[i] = 'X';
			} else if (p->type == 3) {
				k = r.tb > r.qb? r.tb - r.qb : 0;
				for (i = k; i < str.l; ++i) ks->seq.s[i] = 'X';
			}
		}
		putchar(ks->qual.l? '@' : '>');
		puts(ks->name.s);
		puts(ks->seq.s);
		if (ks->qual.l) {
			puts("+");
			puts(ks->qual.s);
		}
	}
	free(str.s);
	kseq_destroy(ks);
	gzclose(fp);

	for (j = 0; j < n_adaps; ++j) {
		ta_adap_t *p = &adaps[j];
		fprintf(stderr, "%-15ld ", (long)p->cnt);
		for (i = 0; i < p->len; ++i) fputc("ACGTN"[(int)p->seq[i]], stderr);
		fputc('\n', stderr);
		free(p->seq);
		free(p->qp);
	}
	free(adaps);
	return 0;
}
