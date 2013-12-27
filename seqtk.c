/* The MIT License

   Copyright (c) 20082-2012 by Heng Li <lh3@me.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int n, m;
	uint64_t *a;
} reglist_t;

#include "khash.h"
KHASH_MAP_INIT_STR(reg, reglist_t)

typedef kh_reg_t reghash_t;

reghash_t *stk_reg_read(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t *str;
	// read the list
	str = calloc(1, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int beg = -1, end = -1;
		reglist_t *p;
		khint_t k = kh_get(reg, h, str->s);
		if (k == kh_end(h)) {
			int ret;
			char *s = strdup(str->s);
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(reglist_t));
		}
		p = &kh_val(h, k);
		if (dret != '\n') {
			if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
				beg = atoi(str->s);
				if (dret != '\n') {
					if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
						end = atoi(str->s);
						if (end < 0) end = -1;
					}
				}
			}
		}
		// skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (end < 0 && beg > 0) end = beg, beg = beg - 1; // if there is only one column
		if (beg < 0) beg = 0, end = INT_MAX;
		if (p->n == p->m) {
			p->m = p->m? p->m<<1 : 4;
			p->a = realloc(p->a, p->m * 8);
		}
		p->a[p->n++] = (uint64_t)beg<<32 | end;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	return h;
}

static void cpy_kstr(kstring_t *dst, const kstring_t *src)
{
	if (src->l == 0) return;
	if (src->l + 1 > dst->m) {
		dst->m = src->l + 1;
		kroundup32(dst->m);
		dst->s = realloc(dst->s, dst->m);
	}
	dst->l = src->l;
	memcpy(dst->s, src->s, src->l + 1);
}

static void cpy_kseq(kseq_t *dst, const kseq_t *src)
{
	cpy_kstr(&dst->name, &src->name);
	cpy_kstr(&dst->seq,  &src->seq);
	cpy_kstr(&dst->qual, &src->qual);
	cpy_kstr(&dst->comment, &src->comment);
}

void stk_reg_destroy(reghash_t *h)
{
	khint_t k;
	if (h == 0) return;
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			free(kh_val(h, k).a);
			free((char*)kh_key(h, k));
		}
	}
	kh_destroy(reg, h);
}

/* constant table */

unsigned char seq_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

char *seq_nt16_rev_table = "XACMGRSVTWYHKDBN";
unsigned char seq_nt16to4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
unsigned char seq_nt16comp_table[] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

static void stk_printstr(const kstring_t *s, unsigned line_len)
{
	if (line_len != UINT_MAX) {
		int i, rest = s->l;
		for (i = 0; i < s->l; i += line_len, rest -= line_len) {
			putchar('\n');
			if (rest > line_len) fwrite(s->s + i, 1, line_len, stdout);
			else fwrite(s->s + i, 1, rest, stdout);
		}
		putchar('\n');
	} else {
		putchar('\n');
		puts(s->s);
	}
}

void stk_printseq(const kseq_t *s, int line_len)
{
	putchar(s->qual.l? '@' : '>');
	fputs(s->name.s, stdout);
	if (s->comment.l) {
		putchar(' '); fputs(s->comment.s, stdout);
	}
	stk_printstr(&s->seq, line_len);
	if (s->qual.l) {
		putchar('+');
		stk_printstr(&s->qual, line_len);
	}
}


/* given an interleaved fastq file, drop sequences without an adjacent pair */
int stk_pairfq(int argc, char *argv[]){
	gzFile fp;
	kseq_t *seq1, *seq2;
    int skipped = 0;

    if (argc != 2){
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk pairfq <interleaved.fq>\n\n");
		fprintf(stderr, "\n");
        return 1;
    }
    fp = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	seq1 = kseq_init(fp);
    kseq_read(seq1);

    seq2 = malloc(sizeof(kseq_t));
    seq2->f = seq1->f;

    while (kseq_read(seq2) >= 0){
        if (strcmp(seq1->name.s, seq2->name.s) == 0){
            stk_printseq(seq1, UINT_MAX);
            stk_printseq(seq2, UINT_MAX);
            kseq_read(seq1);
        } else {
            // skip 1 of these reads
            cpy_kseq(seq1, seq2);
            skipped++;
        }
    }
    fprintf(stderr, "seqtk: dropped %d singletons\n", skipped);
	kseq_destroy(seq1);
    free(seq2->seq.s); free(seq2->qual.s); free(seq2->name.s); free(seq2);
	gzclose(fp);
    return 0;
}

/* convert paired-end files to interleaved */
int stk_interleave(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;

    if (argc != 3){
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk interleavefqs <r1.fq r2.fq>\n\n");
		fprintf(stderr, "\n");
        return 1;
    }

    fp1 = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	seq1 = kseq_init(fp1);

    fp2 = strcmp(argv[2], "-")? gzopen(argv[2], "r") : gzdopen(fileno(stdin), "r");
	seq2 = kseq_init(fp2);

    while (kseq_read(seq1) >= 0){
        if (kseq_read(seq2) < 0){
            fprintf(stderr, "seqtk error: %s has fewer records than %s", argv[2], argv[1]);
            return 1;
        }
        stk_printseq(seq1, UINT_MAX);
        stk_printseq(seq2, UINT_MAX);
    }
    if (kseq_read(seq2) > 0){
        fprintf(stderr, "seqtk error: %s has fewer records than %s\n", argv[1], argv[2]);
        return 1;
    }
	kseq_destroy(seq1);
	kseq_destroy(seq2);
	gzclose(fp1);
	gzclose(fp2);
    return 0;
}


/* quality based trimming with Mott's algorithm */
int stk_trimfq(int argc, char *argv[])
{ // FIXME: when a record with zero length will always be treated as a fasta record
	gzFile fp;
	kseq_t *seq;
	double param = 0.05, q_int2real[128];
	int i, c, min_len = 30, left = 0, right = 0;
	while ((c = getopt(argc, argv, "l:q:b:e:")) >= 0) {
		switch (c) {
			case 'q': param = atof(optarg); break;
			case 'l': min_len = atoi(optarg); break;
			case 'b': left = atoi(optarg); break;
			case 'e': right = atoi(optarg); break;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk trimfq [options] <in.fq>\n\n");
		fprintf(stderr, "Options: -q FLOAT    error rate threshold (disabled by -b/-e) [%.2f]\n", param);
		fprintf(stderr, "         -l INT      maximally trim down to INT bp (disabled by -b/-e) [%d]\n", min_len);
		fprintf(stderr, "         -b INT      trim INT bp from left (non-zero to disable -q/-l) [0]\n");
		fprintf(stderr, "         -e INT      trim INT bp from right (non-zero to disable -q/-l) [0]\n");
		fprintf(stderr, "\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	for (i = 0; i < 128; ++i)
		q_int2real[i] = pow(10., -(i - 33) / 10.);
	while (kseq_read(seq) >= 0) {
		int beg, tmp, end;
		double s, max;
		if (left || right) {
			beg = left; end = seq->seq.l - right;
			if (beg >= end) beg = end = 0;
		} else if (seq->qual.l > min_len) {
			for (i = 0, beg = tmp = 0, end = seq->qual.l, s = max = 0.; i < seq->qual.l; ++i) {
				int q = seq->qual.s[i];
				if (q < 36) q = 36;
				if (q > 127) q = 127;
				s += param - q_int2real[q];
				if (s > max) max = s, beg = tmp, end = i + 1;
				if (s < 0) s = 0, tmp = i + 1;
			}

			/* max never set; all low qual, just give first min_len bp */
			if (max == 0.) beg = 0, end = min_len;

			if (end - beg < min_len) { // window-based 
				int is, imax;
				for (i = 0, is = 0; i < min_len; ++i)
					is += seq->qual.s[i] - 33;
				for (imax = is, beg = 0; i < seq->qual.l; ++i) {
					is += (int)seq->qual.s[i] - seq->qual.s[i - min_len];
					if (imax < is) imax = is, beg = i - min_len + 1;
				}
				end = beg + min_len;
			}
		} else beg = 0, end = seq->seq.l;
		putchar(seq->qual.l? '@' : '>'); fputs(seq->name.s, stdout); 
		if (seq->comment.l) {
			putchar(' '); puts(seq->comment.s);
		} else putchar('\n');
		fwrite(seq->seq.s + beg, 1, end - beg, stdout); putchar('\n');
		if (seq->qual.l) {
			puts("+");
			fwrite(seq->qual.s + beg, 1, end - beg, stdout); putchar('\n');
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

/* composition */
int stk_comp(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l, c, upper_only = 0;
	reghash_t *h = 0;
	reglist_t dummy;
	while ((c = getopt(argc, argv, "ur:")) >= 0) {
		switch (c) {
			case 'u': upper_only = 1; break;
			case 'r': h = stk_reg_read(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage:  seqtk comp [-u] [-r in.bed] <in.fa>\n\n");
		fprintf(stderr, "Output format: chr, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts\n");
		return 1;
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	seq = kseq_init(fp);
	dummy.n= dummy.m = 1; dummy.a = calloc(1, 8);
	while ((l = kseq_read(seq)) >= 0) {
		int i, k;
		reglist_t *p = 0;
		if (h) {
			khint_t k = kh_get(reg, h, seq->name.s);
			if (k != kh_end(h)) p = &kh_val(h, k);
		} else {
			p = &dummy;
			dummy.a[0] = l;
		}
		for (k = 0; p && k < p->n; ++k) {
			int beg = p->a[k]>>32, end = p->a[k]&0xffffffff;
			int la, lb, lc, na, nb, nc, cnt[11];
			if (beg > 0) la = seq->seq.s[beg-1], lb = seq_nt16_table[la], lc = bitcnt_table[lb];
			else la = 'a', lb = -1, lc = 0;
			na = seq->seq.s[beg]; nb = seq_nt16_table[na]; nc = bitcnt_table[nb];
			memset(cnt, 0, 11 * sizeof(int));
			for (i = beg; i < end; ++i) {
				int is_CpG = 0, a, b, c;
				a = na; b = nb; c = nc;
				na = seq->seq.s[i+1]; nb = seq_nt16_table[na]; nc = bitcnt_table[nb];
				if (b == 2 || b == 10) { // C or Y
					if (nb == 4 || nb == 5) is_CpG = 1;
				} else if (b == 4 || b == 5) { // G or R
					if (lb == 2 || lb == 10) is_CpG = 1;
				}
				if (upper_only == 0 || isupper(a)) {
					if (c > 1) ++cnt[c+2];
					if (c == 1) ++cnt[seq_nt16to4_table[b]];
					if (b == 10 || b == 5) ++cnt[9];
					else if (c == 2) {
						++cnt[8];
					}
					if (is_CpG) {
						++cnt[7];
						if (b == 10 || b == 5) ++cnt[10];
					}
				}
				la = a; lb = b; lc = c;
			}
			if (h) printf("%s\t%d\t%d", seq->name.s, beg, end);
			else printf("%s\t%d", seq->name.s, l);
			for (i = 0; i < 11; ++i) printf("\t%d", cnt[i]);
			putchar('\n');
		}
		fflush(stdout);
	}
	free(dummy.a);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int stk_randbase(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: seqtk randbase <in.fa>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		int i;
		printf(">%s", seq->name.s);
		for (i = 0; i < l; ++i) {
			int c, b, a, j, k, m;
			b = seq->seq.s[i];
			c = seq_nt16_table[b];
			a = bitcnt_table[c];
			if (a == 2) {
				m = (drand48() < 0.5);
				for (j = k = 0; j < 4; ++j) {
					if ((1<<j & c) == 0) continue;
					if (k == m) break;
					++k;
				}
				seq->seq.s[i] = islower(b)? "acgt"[j] : "ACGT"[j];
			}
			if (i%60 == 0) putchar('\n');
			putchar(seq->seq.s[i]);
		}
		putchar('\n');
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int stk_hety(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l, c, win_size = 50000, n_start = 5, win_step, is_lower_mask = 0;
	char *buf;
	uint32_t cnt[3];
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk hety [options] <in.fa>\n\n");
		fprintf(stderr, "Options: -w INT   window size [%d]\n", win_size);
		fprintf(stderr, "         -t INT   # start positions in a window [%d]\n", n_start);
		fprintf(stderr, "         -m       treat lowercases as masked\n");
		fprintf(stderr, "\n");
		return 1;
	}
	while ((c = getopt(argc, argv, "w:t:m")) >= 0) {
		switch (c) {
		case 'w': win_size = atoi(optarg); break;
		case 't': n_start = atoi(optarg); break;
		case 'm': is_lower_mask = 1; break;
		}
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	seq = kseq_init(fp);
	win_step = win_size / n_start;
	buf = calloc(win_size, 1);
	while ((l = kseq_read(seq)) >= 0) {
		int x, i, y, z, next = 0;
		cnt[0] = cnt[1] = cnt[2] = 0;
		for (i = 0; i <= l; ++i) {
			if ((i >= win_size && i % win_step == 0) || i == l) {
				if (i == l && l >= win_size) {
					for (y = l - win_size; y < next; ++y) --cnt[(int)buf[y % win_size]];
				}
				if (cnt[1] + cnt[2] > 0)
					printf("%s\t%d\t%d\t%.2lf\t%d\t%d\n", seq->name.s, next, i,
						   (double)cnt[2] / (cnt[1] + cnt[2]) * win_size, cnt[1] + cnt[2], cnt[2]);
				next = i;
			}
			if (i < l) {
				y = i % win_size;
				c = seq->seq.s[i];
				if (is_lower_mask && islower(c)) c = 'N';
				c = seq_nt16_table[c];
				x = bitcnt_table[c];
				if (i >= win_size) --cnt[(int)buf[y]];
				buf[y] = z = x > 2? 0 : x == 2? 2 : 1;
				++cnt[z];
			}
		}
	}
	free(buf);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

/* subseq */

int stk_subseq(int argc, char *argv[])
{
	khash_t(reg) *h = kh_init(reg);
	gzFile fp;
	kseq_t *seq;
	int l, i, j, c, is_tab = 0, line = 1024;
	khint_t k;
	while ((c = getopt(argc, argv, "tl:")) >= 0) {
		switch (c) {
		case 't': is_tab = 1; break;
		case 'l': line = atoi(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk subseq [options] <in.fa> <in.bed>|<name.list>\n\n");
		fprintf(stderr, "Options: -t       TAB delimited output\n");
		fprintf(stderr, "         -l INT   sequence line length [%d]\n\n", line);
		fprintf(stderr, "Note: Use 'samtools faidx' if only a few regions are intended.\n\n");
		return 1;
	}
	h = stk_reg_read(argv[optind+1]);
	// subseq
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		reglist_t *p;
		k = kh_get(reg, h, seq->name.s);
		if (k == kh_end(h)) continue;
		p = &kh_val(h, k);
		for (i = 0; i < p->n; ++i) {
			int beg = p->a[i]>>32, end = p->a[i];
			if (beg >= seq->seq.l) {
				fprintf(stderr, "[subseq] %s: %d >= %ld\n", seq->name.s, beg, seq->seq.l);
				continue;
			}
			if (end > seq->seq.l) end = seq->seq.l;
			if (is_tab == 0) {
				printf("%c%s", seq->qual.l == seq->seq.l? '@' : '>', seq->name.s);
				if (beg > 0 || (int)p->a[i] != INT_MAX) {
					if (end == INT_MAX) {
						if (beg) printf(":%d", beg+1);
					} else printf(":%d-%d", beg+1, end);
				}
			} else printf("%s\t%d\t", seq->name.s, beg + 1);
			if (end > seq->seq.l) end = seq->seq.l;
			for (j = 0; j < end - beg; ++j) {
				if (is_tab == 0 && j % line == 0) putchar('\n');
				putchar(seq->seq.s[j + beg]);
			}
			putchar('\n');
			if (seq->qual.l != seq->seq.l || is_tab) continue;
			printf("+");
			for (j = 0; j < end - beg; ++j) {
				if (j % line == 0) putchar('\n');
				putchar(seq->qual.s[j + beg]);
			}
			putchar('\n');
		}
	}
	// free
	kseq_destroy(seq);
	gzclose(fp);
	stk_reg_destroy(h);
	return 0;
}

/* mergefa */
int stk_mergefa(int argc, char *argv[])
{
	gzFile fp[2];
	kseq_t *seq[2];
	int i, l, c, is_intersect = 0, is_haploid = 0, qual = 0, is_mask = 0, is_randhet = 0;
	uint64_t cnt[5];
	while ((c = getopt(argc, argv, "himrq:")) >= 0) {
		switch (c) {
			case 'i': is_intersect = 1; break;
			case 'h': is_haploid = 1; break;
			case 'm': is_mask = 1; break;
			case 'r': is_randhet = 1; break;
			case 'q': qual = atoi(optarg); break;
		}
	}
	if (is_mask && is_intersect) {
		fprintf(stderr, "[%s] `-i' and `-h' cannot be applied at the same time.\n", __func__);
		return 1;
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\nUsage: seqtk mergefa [options] <in1.fa> <in2.fa>\n\n");
		fprintf(stderr, "Options: -q INT   quality threshold [0]\n");
		fprintf(stderr, "         -i       take intersection\n");
		fprintf(stderr, "         -m       convert to lowercase when one of the input base is N\n");
		fprintf(stderr, "         -r       pick a random allele from het\n");
		fprintf(stderr, "         -h       suppress hets in the input\n\n");
		return 1;
	}
	for (i = 0; i < 2; ++i) {
		fp[i] = strcmp(argv[optind+i], "-")? gzopen(argv[optind+i], "r") : gzdopen(fileno(stdin), "r");
		seq[i] = kseq_init(fp[i]);
	}
	cnt[0] = cnt[1] = cnt[2] = cnt[3] = cnt[4] = 0;
	srand48(11);
	while (kseq_read(seq[0]) >= 0) {
		int min_l, c[2], b[2], is_upper;
		kseq_read(seq[1]);
		if (strcmp(seq[0]->name.s, seq[1]->name.s))
			fprintf(stderr, "[%s] Different sequence names: %s != %s\n", __func__, seq[0]->name.s, seq[1]->name.s);
		if (seq[0]->seq.l != seq[1]->seq.l)
			fprintf(stderr, "[%s] Unequal sequence length: %ld != %ld\n", __func__, seq[0]->seq.l, seq[1]->seq.l);
		min_l = seq[0]->seq.l < seq[1]->seq.l? seq[0]->seq.l : seq[1]->seq.l;
		printf(">%s", seq[0]->name.s);
		for (l = 0; l < min_l; ++l) {
			c[0] = seq[0]->seq.s[l]; c[1] = seq[1]->seq.s[l];
			if (seq[0]->qual.l && seq[0]->qual.s[l] - 33 < qual) c[0] = tolower(c[0]);
			if (seq[1]->qual.l && seq[1]->qual.s[l] - 33 < qual) c[1] = tolower(c[1]);
			if (is_intersect) is_upper = (isupper(c[0]) || isupper(c[1]))? 1 : 0;
			else if (is_mask) is_upper = (isupper(c[0]) || isupper(c[1]))? 1 : 0;
			else is_upper = (isupper(c[0]) && isupper(c[1]))? 1 : 0;
			c[0] = seq_nt16_table[c[0]]; c[1] = seq_nt16_table[c[1]];
			if (c[0] == 0) c[0] = 15;
			if (c[1] == 0) c[1] = 15;
			b[0] = bitcnt_table[c[0]];
			b[1] = bitcnt_table[c[1]];
			if (is_upper) {
				if (b[0] == 1 && b[1] == 1) {
					if (c[0] == c[1]) ++cnt[0];
					else ++cnt[1];
				} else if (b[0] == 1 && b[1] == 2) ++cnt[2];
				else if (b[0] == 2 && b[1] == 1) ++cnt[3];
				else if (b[0] == 2 && b[1] == 2) ++cnt[4];
			}
			if (is_haploid && (b[0] > 1 || b[1] > 1)) is_upper = 0;
			if (is_intersect) {
				c[0] = c[0] & c[1];
				if (c[0] == 0) is_upper = 0; // FIXME: is this a bug - c[0] cannot be 0!
			} else if (is_mask) {
				if (c[0] == 15 || c[1] == 15) is_upper = 0;
				c[0] &= c[1];
				if (c[0] == 0) is_upper = 0;
			} else if (is_randhet) {
				if (b[0] == 1 && b[1] == 1) { // two homs
					c[0] |= c[1];
				} else if (((b[0] == 1 && b[1] == 2) || (b[0] == 2 && b[1] == 1)) && (c[0]&c[1])) { // one hom, one het
					c[0] = (lrand48()&1)? (c[0] & c[1]) : (c[0] | c[1]);
				} else if (b[0] == 2 && b[1] == 2 && c[0] == c[1]) { // double hets
					if (lrand48()&1) {
						if (lrand48()&1) {
							for (i = 8; i >= 1; i >>= 1) // pick the "larger" allele
								if (c[0]&i) c[0] &= i;
						} else {
							for (i = 1; i <= 8; i <<= 1) // pick the "smaller" allele
								if (c[0]&i) c[0] &= i;
						}
					} // else set as het
				} else is_upper = 0;
			} else c[0] |= c[1];
			c[0] = seq_nt16_rev_table[c[0]];
			if (!is_upper) c[0] = tolower(c[0]);
			if (l%60 == 0) putchar('\n');
			putchar(c[0]);
		}
		putchar('\n');
	}
	fprintf(stderr, "[%s] (same,diff,hom-het,het-hom,het-het)=(%ld,%ld,%ld,%ld,%ld)\n", __func__, (long)cnt[0], (long)cnt[1], (long)cnt[2], (long)cnt[3], (long)cnt[4]);
	return 0;
}

int stk_famask(int argc, char *argv[])
{
	gzFile fp[2];
	kseq_t *seq[2];
	int i, l;
	if (argc < 3) {
		fprintf(stderr, "Usage: seqtk famask <src.fa> <mask.fa>\n");
		return 1;
	}
	for (i = 0; i < 2; ++i) {
		fp[i] = strcmp(argv[optind+i], "-")? gzopen(argv[optind+i], "r") : gzdopen(fileno(stdin), "r");
		seq[i] = kseq_init(fp[i]);
	}
	while (kseq_read(seq[0]) >= 0) {
		int min_l, c[2];
		kseq_read(seq[1]);
		if (strcmp(seq[0]->name.s, seq[1]->name.s))
			fprintf(stderr, "[%s] Different sequence names: %s != %s\n", __func__, seq[0]->name.s, seq[1]->name.s);
		if (seq[0]->seq.l != seq[1]->seq.l)
			fprintf(stderr, "[%s] Unequal sequence length: %ld != %ld\n", __func__, seq[0]->seq.l, seq[1]->seq.l);
		min_l = seq[0]->seq.l < seq[1]->seq.l? seq[0]->seq.l : seq[1]->seq.l;
		printf(">%s", seq[0]->name.s);
		for (l = 0; l < min_l; ++l) {
			c[0] = seq[0]->seq.s[l]; c[1] = seq[1]->seq.s[l];
			if (c[1] == 'x') c[0] = tolower(c[0]);
			else if (c[1] != 'X') c[0] = c[1];
			if (l%60 == 0) putchar('\n');
			putchar(c[0]);
		}
		putchar('\n');
	}
	return 0;
}

int stk_mutfa(int argc, char *argv[])
{
	khash_t(reg) *h = kh_init(reg);
	gzFile fp;
	kseq_t *seq;
	kstream_t *ks;
	int l, i, dret;
	kstring_t *str;
	khint_t k;
	if (argc < 3) {
		fprintf(stderr, "Usage: seqtk mutfa <in.fa> <in.snp>\n\n");
		fprintf(stderr, "Note: <in.snp> contains at least four columns per line which are:\n");
		fprintf(stderr, "      'chr  1-based-pos  any  base-changed-to'.\n");
		return 1;
	}
	// read the list
	str = calloc(1, sizeof(kstring_t));
	fp = strcmp(argv[2], "-")? gzopen(argv[2], "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		char *s = strdup(str->s);
		int beg = 0, ret;
		reglist_t *p;
		k = kh_get(reg, h, s);
		if (k == kh_end(h)) {
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(reglist_t));
		}
		p = &kh_val(h, k);
		if (ks_getuntil(ks, 0, str, &dret) > 0) beg = atol(str->s) - 1; // 2nd col
		ks_getuntil(ks, 0, str, &dret); // 3rd col
		ks_getuntil(ks, 0, str, &dret); // 4th col
		// skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (isalpha(str->s[0]) && str->l == 1) {
			if (p->n == p->m) {
				p->m = p->m? p->m<<1 : 4;
				p->a = realloc(p->a, p->m * 8);
			}
			p->a[p->n++] = (uint64_t)beg<<32 | str->s[0];
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	// mutfa
	fp = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		reglist_t *p;
		k = kh_get(reg, h, seq->name.s);
		if (k != kh_end(h)) {
			p = &kh_val(h, k);
			for (i = 0; i < p->n; ++i) {
				int beg = p->a[i]>>32;
				if (beg < seq->seq.l)
					seq->seq.s[beg] = (int)p->a[i];
			}
		}
		printf(">%s", seq->name.s);
		for (i = 0; i < l; ++i) {
			if (i%60 == 0) putchar('\n');
			putchar(seq->seq.s[i]);
		}
		putchar('\n');
	}
	// free
	kseq_destroy(seq);
	gzclose(fp);
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			free(kh_val(h, k).a);
			free((char*)kh_key(h, k));
		}
	}
	kh_destroy(reg, h);
	return 0;
}

int stk_listhet(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int i, l;
	if (argc == 1) {
		fprintf(stderr, "Usage: seqtk listhet <in.fa>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		for (i = 0; i < l; ++i) {
			int b = seq->seq.s[i];
			if (bitcnt_table[seq_nt16_table[b]] == 2)
				printf("%s\t%d\t%c\n", seq->name.s, i+1, b);
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

/* cutN */
static int cutN_min_N_tract = 1000;
static int cutN_nonN_penalty = 10;

static int find_next_cut(const kseq_t *ks, int k, int *begin, int *end)
{
	int i, b, e;
	while (k < ks->seq.l) {
		if (seq_nt16_table[(int)ks->seq.s[k]] == 15) {
			int score, max;
			score = 0; e = max = -1;
			for (i = k; i < ks->seq.l && score >= 0; ++i) { /* forward */
				if (seq_nt16_table[(int)ks->seq.s[i]] == 15) ++score;
				else score -= cutN_nonN_penalty;
				if (score > max) max = score, e = i;
			}
			score = 0; b = max = -1;
			for (i = e; i >= 0 && score >= 0; --i) { /* backward */
				if (seq_nt16_table[(int)ks->seq.s[i]] == 15) ++score;
				else score -= cutN_nonN_penalty;
				if (score > max) max = score, b = i;
			}
			if (e + 1 - b >= cutN_min_N_tract) {
				*begin = b;
				*end = e + 1;
				return *end;
			}
			k = e + 1;
		} else ++k;
	}
	return -1;
}
static void print_seq(FILE *fpout, const kseq_t *ks, int begin, int end)
{
	int i;
	if (begin >= end) return; // FIXME: why may this happen? Understand it!
	fprintf(fpout, ">%s:%d-%d", ks->name.s, begin+1, end);
	for (i = begin; i < end && i < ks->seq.l; ++i) {
		if ((i - begin)%60 == 0) fputc('\n', fpout);
		fputc(ks->seq.s[i], fpout);
	}
	fputc('\n', fpout);
}
int stk_cutN(int argc, char *argv[])
{
	int c, l, gap_only = 0;
	gzFile fp;
	kseq_t *ks;
	while ((c = getopt(argc, argv, "n:p:g")) >= 0) {
		switch (c) {
		case 'n': cutN_min_N_tract = atoi(optarg); break;
		case 'p': cutN_nonN_penalty = atoi(optarg); break;
		case 'g': gap_only = 1; break;
		default: return 1;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk cutN [options] <in.fa>\n\n");
		fprintf(stderr, "Options: -n INT    min size of N tract [%d]\n", cutN_min_N_tract);
		fprintf(stderr, "         -p INT    penalty for a non-N [%d]\n", cutN_nonN_penalty);
		fprintf(stderr, "         -g        print gaps only, no sequence\n\n");
		return 1;
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	ks = kseq_init(fp);
	while ((l = kseq_read(ks)) >= 0) {
		int k = 0, begin = 0, end = 0;
		while (find_next_cut(ks, k, &begin, &end) >= 0) {
			if (begin != 0) {
				if (gap_only) printf("%s\t%d\t%d\n", ks->name.s, begin, end);
				else print_seq(stdout, ks, k, begin);
			}
			k = end;
		}
		if (!gap_only) print_seq(stdout, ks, k, l);
	}
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}

int stk_hrun(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *ks;
	int min_len = 7, l = 0, c = 0, beg = 0, i;
	if (argc == optind) {
		fprintf(stderr, "Usage: seqtk hrun <in.fa> [minLen=%d]\n", min_len);
		return 1;
	}
	if (argc == optind + 2) min_len = atoi(argv[optind+1]);
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		c = ks->seq.s[0]; l = 1; beg = 0;
		for (i = 1; i < ks->seq.l; ++i) {
			if (ks->seq.s[i] != c) {
				if (l >= min_len) printf("%s\t%d\t%d\t%c\n", ks->name.s, beg, beg + l, c);
				c = ks->seq.s[i]; l = 1; beg = i;
			} else ++l;
		}
	}
	if (l >= min_len) printf("%s\t%d\t%d\t%c\n", ks->name.s, beg, beg + l, c);
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}

/* sample */

int stk_sample(int argc, char *argv[])
{
	int c;
	uint64_t i, num = 0, n_seqs = 0;
	double frac = 0.;
	gzFile fp;
	kseq_t *seq, *buf = 0;

	srand48(11);
	while ((c = getopt(argc, argv, "s:")) >= 0)
		switch (c) {
		case 's': srand48(atoi(optarg)); break;
		}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: seqtk sample [-s seed=11] <in.fa> <frac>|<number>\n\n");
		fprintf(stderr, "Warning: Large memory consumption for large <number>.\n");
		return 1;
	}
	frac = atof(argv[optind+1]);
	if (frac > 1.) num = (uint64_t)(frac + .499), frac = 0.;
	if (num > 0) buf = calloc(num, sizeof(kseq_t));

	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		double r = drand48();
		++n_seqs;
		if (num) {
			uint64_t y = n_seqs - 1 < num? n_seqs - 1 : (uint64_t)(r * n_seqs);
			if (y < num) cpy_kseq(&buf[y], seq);
		} else if (r < frac) stk_printseq(seq, UINT_MAX);
	}
	kseq_destroy(seq);
	gzclose(fp);
	for (i = 0; i < num; ++i) {
		kseq_t *p = &buf[i];
		if (p->seq.l) stk_printseq(p, UINT_MAX);
		free(p->seq.s); free(p->qual.s); free(p->name.s);
	}
	free(buf);
	return 0;
}

/* seq */

void stk_mask(kseq_t *seq, const khash_t(reg) *h, int is_complement, int mask_chr)
{
	unsigned i, j;
	khiter_t k;
	k = kh_get(reg, h, seq->name.s);
	if (k == kh_end(h)) { // not found in the hash table
		if (is_complement) {
			if (mask_chr) {
				for (j = 0; j < seq->seq.l; ++j)
					seq->seq.s[j] = mask_chr;
			} else {
				for (j = 0; j < seq->seq.l; ++j)
					seq->seq.s[j] = tolower(seq->seq.s[j]);
			}
		}
	} else {
		reglist_t *p = &kh_val(h, k);
		if (!is_complement) {
			for (i = 0; i < p->n; ++i) {
				unsigned beg = p->a[i]>>32, end = p->a[i];
				if (beg >= seq->seq.l) continue;
				if (end > seq->seq.l) end = seq->seq.l;
				if (!mask_chr) for (j = beg; j < end; ++j) seq->seq.s[j] = tolower(seq->seq.s[j]);
				else for (j = beg; j < end; ++j) seq->seq.s[j] = mask_chr;
			}
		} else {
			int8_t *mask = calloc(seq->seq.l, 1);
			for (i = 0; i < p->n; ++i) {
				unsigned beg = p->a[i]>>32, end = p->a[i];
				if (end >= seq->seq.l) end = seq->seq.l;
				for (j = beg; j < end; ++j) mask[j] = 1;
			}
			if (mask_chr) {
				for (j = 0; j < seq->seq.l; ++j)
					if (mask[j] == 0) seq->seq.s[j] = mask_chr;
			} else {
				for (j = 0; j < seq->seq.l; ++j)
					if (mask[j] == 0) seq->seq.s[j] = tolower(seq->seq.s[j]);
			}
			free(mask);
		}
	}
}

int stk_seq(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int c, qual_thres = 0, flag = 0, qual_shift = 33, mask_chr = 0, min_len = 0;
	unsigned i, line_len = 0;
	int64_t n_seqs = 0;
	double frac = 1.;
	khash_t(reg) *h = 0;

	srand48(11);
	while ((c = getopt(argc, argv, "12q:l:Q:aACrn:s:f:M:L:cV")) >= 0) {
		switch (c) {
			case 'a':
			case 'A': flag |= 1; break;
			case 'C': flag |= 2; break;
			case 'r': flag |= 4; break;
			case 'c': flag |= 8; break;
			case '1': flag |= 16; break;
			case '2': flag |= 32; break;
			case 'V': flag |= 64; break;
			case 'M': h = stk_reg_read(optarg); break;
			case 'n': mask_chr = *optarg; break;
			case 'Q': qual_shift = atoi(optarg); break;
			case 'q': qual_thres = atoi(optarg); break;
			case 'l': line_len = atoi(optarg); break;
			case 'L': min_len = atoi(optarg); break;
			case 's': srand48(atoi(optarg)); break;
			case 'f': frac = atof(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk seq [options] <in.fq>|<in.fa>\n\n");
		fprintf(stderr, "Options: -q INT    mask bases with quality lower than INT [0]\n");
		fprintf(stderr, "         -n CHAR   masked bases converted to CHAR; 0 for lowercase [0]\n");
		fprintf(stderr, "         -l INT    number of residues per line; 0 for 2^32-1 [%d]\n", line_len);
		fprintf(stderr, "         -Q INT    quality shift: ASCII-INT gives base quality [%d]\n", qual_shift);
		fprintf(stderr, "         -s INT    random seed (effective with -f) [11]\n");
		fprintf(stderr, "         -f FLOAT  sample FLOAT fraction of sequences [1]\n");
		fprintf(stderr, "         -M FILE   mask regions in BED or name list FILE [null]\n");
		fprintf(stderr, "         -L INT    drop sequences with length shorter than INT [0]\n");
		fprintf(stderr, "         -c        mask complement region (effective with -M)\n");
		fprintf(stderr, "         -r        reverse complement\n");
		fprintf(stderr, "         -A        force FASTA output (discard quality)\n");
		fprintf(stderr, "         -C        drop comments at the header lines\n");
		fprintf(stderr, "         -1        output the 2n-1 reads only\n");
		fprintf(stderr, "         -2        output the 2n reads only\n");
		fprintf(stderr, "         -V        shift quality by '(-Q) - 33'\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (line_len == 0) line_len = UINT_MAX;
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	qual_thres += qual_shift;
	while (kseq_read(seq) >= 0) {
		++n_seqs;
		if (seq->seq.l < min_len) continue; // NB: length filter before taking random
		if (frac < 1. && drand48() >= frac) continue;
		if (flag & 48) { // then choose odd/even reads only
			if ((flag&16) && (n_seqs&1) == 0) continue;
			if ((flag&32) && (n_seqs&1) == 1) continue;
		}
		if (seq->qual.l && qual_thres > qual_shift) {
			if (mask_chr) {
				for (i = 0; i < seq->seq.l; ++i)
					if (seq->qual.s[i] < qual_thres)
						seq->seq.s[i] = mask_chr;
			} else {
				for (i = 0; i < seq->seq.l; ++i)
					if (seq->qual.s[i] < qual_thres)
						seq->seq.s[i] = tolower(seq->seq.s[i]);
			}
		}
		if (flag & 1) seq->qual.l = 0;
		if (flag & 2) seq->comment.l = 0;
		if (h) stk_mask(seq, h, flag&8, mask_chr); // masking
		if (flag & 4) { // reverse complement
			int c0, c1;
			for (i = 0; i < seq->seq.l>>1; ++i) { // reverse complement sequence
				c0 = comp_tab[(int)seq->seq.s[i]];
				c1 = comp_tab[(int)seq->seq.s[seq->seq.l - 1 - i]];
				seq->seq.s[i] = c1;
				seq->seq.s[seq->seq.l - 1 - i] = c0;
			}
			if (seq->seq.l & 1) // complement the remaining base
				seq->seq.s[seq->seq.l>>1] = comp_tab[(int)seq->seq.s[seq->seq.l>>1]];
			if (seq->qual.l) {
				for (i = 0; i < seq->seq.l>>1; ++i) // reverse quality
					c0 = seq->qual.s[i], seq->qual.s[i] = seq->qual.s[seq->qual.l - 1 - i], seq->qual.s[seq->qual.l - 1 - i] = c0;
			}
		}
		if ((flag & 64) && seq->qual.l && qual_shift != 33)
			for (i = 0; i < seq->qual.l; ++i)
				seq->qual.s[i] -= qual_shift - 33;
		stk_printseq(seq, line_len);
	}
	kseq_destroy(seq);
	gzclose(fp);
	stk_reg_destroy(h);
	return 0;
}

/* main function */
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   seqtk <command> <arguments>\n");
	fprintf(stderr, "Version: 1.0-r32\n\n");
	fprintf(stderr, "Command: seq        common transformation of FASTA/Q\n");
	fprintf(stderr, "         comp       get the nucleotide composition of FASTA/Q\n");
	fprintf(stderr, "         sample     subsample sequences\n");
	fprintf(stderr, "         interleave interleave paired-end fastq sequences\n");
	fprintf(stderr, "         pairfq     drop un-paired reads from interleaved FASTQ\n");
	fprintf(stderr, "         subseq     extract subsequences from FASTA/Q\n");
	fprintf(stderr, "         trimfq     trim FASTQ using the Phred algorithm\n\n");
	fprintf(stderr, "         hety       regional heterozygosity\n");
	fprintf(stderr, "         mutfa      point mutate FASTA at specified positions\n");
	fprintf(stderr, "         mergefa    merge two FASTA/Q files\n");
	fprintf(stderr, "         randbase   choose a random base from hets\n");
	fprintf(stderr, "         cutN       cut sequence at long N\n");
	fprintf(stderr, "         listhet    extract the position of each het\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();
	if (strcmp(argv[1], "comp") == 0) stk_comp(argc-1, argv+1);
	else if (strcmp(argv[1], "hety") == 0) stk_hety(argc-1, argv+1);
	else if (strcmp(argv[1], "subseq") == 0) stk_subseq(argc-1, argv+1);
	else if (strcmp(argv[1], "interleave") == 0) stk_interleave(argc-1, argv+1);
	else if (strcmp(argv[1], "pairfq") == 0) stk_pairfq(argc-1, argv+1);
	else if (strcmp(argv[1], "mutfa") == 0) stk_mutfa(argc-1, argv+1);
	else if (strcmp(argv[1], "mergefa") == 0) stk_mergefa(argc-1, argv+1);
	else if (strcmp(argv[1], "randbase") == 0) stk_randbase(argc-1, argv+1);
	else if (strcmp(argv[1], "cutN") == 0) stk_cutN(argc-1, argv+1);
	else if (strcmp(argv[1], "listhet") == 0) stk_listhet(argc-1, argv+1);
	else if (strcmp(argv[1], "famask") == 0) stk_famask(argc-1, argv+1);
	else if (strcmp(argv[1], "trimfq") == 0) stk_trimfq(argc-1, argv+1);
	else if (strcmp(argv[1], "hrun") == 0) stk_hrun(argc-1, argv+1);
	else if (strcmp(argv[1], "sample") == 0) stk_sample(argc-1, argv+1);
	else if (strcmp(argv[1], "seq") == 0) stk_seq(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized commad '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
