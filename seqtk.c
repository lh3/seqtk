/* The MIT License

   Copyright (c) 2008-2016 Broad Institute

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
#include <inttypes.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int n, m;
	uint64_t *a;
} reglist_t;

#include "khash.h"
KHASH_MAP_INIT_STR(reg, reglist_t)
KHASH_SET_INIT_INT64(64)

typedef kh_reg_t reghash_t;

reghash_t *stk_reg_read(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t *str;
	// read the list
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	str = calloc(1, sizeof(kstring_t));
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

unsigned char seq_nt6_table[256] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
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
	if (line_len != UINT_MAX && line_len != 0) {
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

static inline void stk_printseq_renamed(const kseq_t *s, int line_len, const char *prefix, int64_t n)
{
	putchar(s->qual.l? '@' : '>');
	if (n >= 0) {
		if (prefix) fputs(prefix, stdout);
		printf("%lld", (long long)n);
	} else fputs(s->name.s, stdout);
	if (s->comment.l) {
		putchar(' '); fputs(s->comment.s, stdout);
	}
	stk_printstr(&s->seq, line_len);
	if (s->qual.l) {
		putchar('+');
		stk_printstr(&s->qual, line_len);
	}
}

static inline void stk_printseq(const kseq_t *s, int line_len)
{
	stk_printseq_renamed(s, line_len, 0, -1);
}

/* 
   64-bit Mersenne Twister pseudorandom number generator. Adapted from:

     http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c

   which was written by Takuji Nishimura and Makoto Matsumoto and released
   under the 3-clause BSD license.
*/

typedef uint64_t krint64_t;

struct _krand_t;
typedef struct _krand_t krand_t;

#define KR_NN 312
#define KR_MM 156
#define KR_UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define KR_LM 0x7FFFFFFFULL /* Least significant 31 bits */

struct _krand_t {
	int mti;
	krint64_t mt[KR_NN];
};

static void kr_srand0(krint64_t seed, krand_t *kr)
{
	kr->mt[0] = seed;
	for (kr->mti = 1; kr->mti < KR_NN; ++kr->mti) 
		kr->mt[kr->mti] = 6364136223846793005ULL * (kr->mt[kr->mti - 1] ^ (kr->mt[kr->mti - 1] >> 62)) + kr->mti;
}

krand_t *kr_srand(krint64_t seed)
{
	krand_t *kr;
	kr = malloc(sizeof(krand_t));
	kr_srand0(seed, kr);
	return kr;
}

krint64_t kr_rand(krand_t *kr)
{
	krint64_t x;
	static const krint64_t mag01[2] = { 0, 0xB5026F5AA96619E9ULL };
    if (kr->mti >= KR_NN) {
		int i;
		if (kr->mti == KR_NN + 1) kr_srand0(5489ULL, kr);
        for (i = 0; i < KR_NN - KR_MM; ++i) {
            x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
            kr->mt[i] = kr->mt[i + KR_MM] ^ (x>>1) ^ mag01[(int)(x&1)];
        }
        for (; i < KR_NN - 1; ++i) {
            x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
            kr->mt[i] = kr->mt[i + (KR_MM - KR_NN)] ^ (x>>1) ^ mag01[(int)(x&1)];
        }
        x = (kr->mt[KR_NN - 1] & KR_UM) | (kr->mt[0] & KR_LM);
        kr->mt[KR_NN - 1] = kr->mt[KR_MM - 1] ^ (x>>1) ^ mag01[(int)(x&1)];
        kr->mti = 0;
    }
    x = kr->mt[kr->mti++];
    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);
    return x;
}

#define kr_drand(_kr) ((kr_rand(_kr) >> 11) * (1.0/9007199254740992.0))


/* quality based trimming with Mott's algorithm */
int stk_trimfq(int argc, char *argv[])
{ // FIXME: when a record with zero length will always be treated as a fasta record
	gzFile fp;
	kseq_t *seq;
	double param = 0.05, q_int2real[128];
	int i, c, min_len = 30, left = 0, right = 0, fixed_len = -1;
	while ((c = getopt(argc, argv, "l:q:b:e:L:")) >= 0) {
		switch (c) {
			case 'q': param = atof(optarg); break;
			case 'l': min_len = atoi(optarg); break;
			case 'b': left = atoi(optarg); break;
			case 'e': right = atoi(optarg); break;
			case 'L': fixed_len = atoi(optarg); break;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk trimfq [options] <in.fq>\n\n");
		fprintf(stderr, "Options: -q FLOAT    error rate threshold (disabled by -b/-e) [%.2f]\n", param);
		fprintf(stderr, "         -l INT      maximally trim down to INT bp (disabled by -b/-e) [%d]\n", min_len);
		fprintf(stderr, "         -b INT      trim INT bp from left (non-zero to disable -q/-l) [0]\n");
		fprintf(stderr, "         -e INT      trim INT bp from right (non-zero to disable -q/-l) [0]\n");
		fprintf(stderr, "         -L INT      retain at most INT bp from the 5'-end (non-zero to disable -q/-l) [0]\n");
		fprintf(stderr, "         -Q          force FASTQ output\n");
		fprintf(stderr, "\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	seq = kseq_init(fp);
	for (i = 0; i < 128; ++i)
		q_int2real[i] = pow(10., -(i - 33) / 10.);
	while (kseq_read(seq) >= 0) {
		int beg, tmp, end;
		double s, max;
		if (left || right || fixed_len > 0) {
			beg = left; end = seq->seq.l - right;
			if (beg >= end) beg = end = 0;
			if (fixed_len > 0 && end - beg > fixed_len) end = beg + fixed_len;
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
		putchar(seq->is_fastq? '@' : '>'); fputs(seq->name.s, stdout); 
		if (seq->comment.l) {
			putchar(' '); puts(seq->comment.s);
		} else putchar('\n');
		fwrite(seq->seq.s + beg, 1, end - beg, stdout); putchar('\n');
		if (seq->is_fastq) {
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
	if (argc == optind && isatty(fileno(stdin))) {
		fprintf(stderr, "Usage:  seqtk comp [-u] [-r in.bed] <in.fa>\n\n");
		fprintf(stderr, "Output format: chr, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts\n");
		return 1;
	}
	fp = optind < argc && strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
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
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
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
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
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
	int l, i, j, c, is_tab = 0, line = 0;
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
	if (h == 0) {
		fprintf(stderr, "[E::%s] failed to read the list of regions in file '%s'\n", __func__, argv[optind+1]);
		return 1;
	}
	// subseq
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream\n", __func__);
		return 1;
	}
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
				if (seq->comment.l) printf(" %s", seq->comment.s);
			} else printf("%s\t%d\t", seq->name.s, beg + 1);
			if (end > seq->seq.l) end = seq->seq.l;
			for (j = 0; j < end - beg; ++j) {
				if (is_tab == 0 && (j == 0 || (line > 0 && j % line == 0))) putchar('\n');
				putchar(seq->seq.s[j + beg]);
			}
			putchar('\n');
			if (seq->qual.l != seq->seq.l || is_tab) continue;
			printf("+");
			for (j = 0; j < end - beg; ++j) {
				if (j == 0 || (line > 0 && j % line == 0)) putchar('\n');
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
	if (fp[0] == 0 || fp[1] == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
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
	int i, l, c;
	while ((c = getopt(argc, argv, "")) >= 0);
	if (argc - optind < 2) {
		fprintf(stderr, "Usage: seqtk famask <src.fa> <mask.fa>\n");
		return 1;
	}
	for (i = 0; i < 2; ++i) {
		fp[i] = strcmp(argv[optind+i], "-")? gzopen(argv[optind+i], "r") : gzdopen(fileno(stdin), "r");
		seq[i] = kseq_init(fp[i]);
	}
	if (fp[0] == 0 || fp[1] == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
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
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
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
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
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
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
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
	fprintf(fpout, "%c%s:%d-%d", ks->qual.l? '@' : '>', ks->name.s, begin+1, end);
	for (i = begin; i < end && i < ks->seq.l; ++i) {
		if ((i - begin)%60 == 0) fputc('\n', fpout);
		fputc(ks->seq.s[i], fpout);
	}
	fputc('\n', fpout);
	if (ks->qual.l == 0) return;
	fputs("+\n", fpout);
	for (i = begin; i < end && i < ks->qual.l; ++i) {
		if ((i - begin)%60 == 0) fputc('\n', fpout);
		fputc(ks->qual.s[i], fpout);
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
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
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
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
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

int stk_sample(int argc, char *argv[])
{
	int c, twopass = 0;
	uint64_t i, num = 0, n_seqs = 0;
	double frac = 0.;
	gzFile fp;
	kseq_t *seq;
	krand_t *kr = 0;

	while ((c = getopt(argc, argv, "2s:")) >= 0)
		if (c == 's') kr = kr_srand(atol(optarg));
		else if (c == '2') twopass = 1;

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk sample [-2] [-s seed=11] <in.fa> <frac>|<number>\n\n");
		fprintf(stderr, "Options: -s INT       RNG seed [11]\n");
		fprintf(stderr, "         -2           2-pass mode: twice as slow but with much reduced memory\n\n");
		return 1;
	}
	frac = atof(argv[optind+1]);
	if (frac > 1.) num = (uint64_t)(frac + .499), frac = 0.;
	else if (twopass) {
		fprintf(stderr, "[W::%s] when sampling a fraction, option -2 is ignored.", __func__);
		twopass = 0;
	}
	if (kr == 0) kr = kr_srand(11);

	if (!twopass) { // the streaming version
		kseq_t *buf = 0;
		if (num > 0) buf = calloc(num, sizeof(kseq_t));
		if (num > 0 && buf == NULL) {
			fprintf(stderr, "[E::%s] Could not allocate enough memory for %" PRIu64 " sequences. Exiting...\n", __func__, num);
			free(kr);
			return 1;
		}

		fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
		if (fp == 0) {
			fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
			return 1;
		}
		seq = kseq_init(fp);
		n_seqs = 0;
		while (kseq_read(seq) >= 0) {
			double r = kr_drand(kr);
			++n_seqs;
			if (num) {
				uint64_t y = n_seqs - 1 < num? n_seqs - 1 : (uint64_t)(r * n_seqs);
				if (y < num) cpy_kseq(&buf[y], seq);
			} else if (r < frac) stk_printseq(seq, UINT_MAX);
		}
		for (i = 0; i < num; ++i) {
			kseq_t *p = &buf[i];
			if (p->seq.l) stk_printseq(p, UINT_MAX);
			free(p->seq.s); free(p->qual.s); free(p->name.s);
		}
		if (buf != NULL) free(buf);
	} else {
		uint64_t *buf;
		khash_t(64) *hash;
		int absent;

		if (strcmp(argv[optind], "-") == 0) {
			fprintf(stderr, "[E::%s] in the 2-pass mode, the input cannot be STDIN.\n", __func__);
			free(kr);
			return 1;
		}

		// 1st pass
		buf = malloc(num * 8);
		for (i = 0; i < num; ++i) buf[i] = UINT64_MAX;
		fp = gzopen(argv[optind], "r");
		if (fp == 0) {
			fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
			return 1;
		}
		seq = kseq_init(fp);
		n_seqs = 0;
		while (kseq_read(seq) >= 0) {
			double r = kr_drand(kr);
			uint64_t y;
			++n_seqs;
			y = n_seqs - 1 < num? n_seqs - 1 : (uint64_t)(r * n_seqs);
			if (y < num) buf[y] = n_seqs;
		}
		kseq_destroy(seq);
		gzclose(fp);
		hash = kh_init(64);
		for (i = 0; i < num; ++i) kh_put(64, hash, buf[i], &absent);
		free(buf);
		// 2nd pass
		fp = gzopen(argv[optind], "r");
		seq = kseq_init(fp);
		n_seqs = 0;
		while (kseq_read(seq) >= 0)
			if (kh_get(64, hash, ++n_seqs) != kh_end(hash))
				stk_printseq(seq, UINT_MAX);
		kh_destroy(64, hash);
	}

	kseq_destroy(seq);
	gzclose(fp);
	free(kr);
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
	int c, qual_thres = 0, flag = 0, qual_shift = 33, mask_chr = 0, min_len = 0, max_q = 255;
	unsigned i, line_len = 0;
	int64_t n_seqs = 0;
	double frac = 1.;
	khash_t(reg) *h = 0;
	krand_t *kr = 0;

	while ((c = getopt(argc, argv, "N12q:l:Q:aACrn:s:f:M:L:cVUX:S")) >= 0) {
		switch (c) {
			case 'a':
			case 'A': flag |= 1; break;
			case 'C': flag |= 2; break;
			case 'r': flag |= 4; break;
			case 'c': flag |= 8; break;
			case '1': flag |= 16; break;
			case '2': flag |= 32; break;
			case 'V': flag |= 64; break;
			case 'N': flag |= 128; break;
			case 'U': flag |= 256; break;
			case 'S': flag |= 512; break;
			case 'M': h = stk_reg_read(optarg); break;
			case 'n': mask_chr = *optarg; break;
			case 'Q': qual_shift = atoi(optarg); break;
			case 'q': qual_thres = atoi(optarg); break;
			case 'X': max_q = atoi(optarg); break;
			case 'l': line_len = atoi(optarg); break;
			case 'L': min_len = atoi(optarg); break;
			case 's': kr = kr_srand(atol(optarg)); break;
			case 'f': frac = atof(optarg); break;
		}
	}
	if (kr == 0) kr = kr_srand(11);
	if (argc == optind && isatty(fileno(stdin))) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk seq [options] <in.fq>|<in.fa>\n\n");
		fprintf(stderr, "Options: -q INT    mask bases with quality lower than INT [0]\n");
		fprintf(stderr, "         -X INT    mask bases with quality higher than INT [255]\n");
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
		fprintf(stderr, "         -N        drop sequences containing ambiguous bases\n");
		fprintf(stderr, "         -1        output the 2n-1 reads only\n");
		fprintf(stderr, "         -2        output the 2n reads only\n");
		fprintf(stderr, "         -V        shift quality by '(-Q) - 33'\n");
		fprintf(stderr, "         -U        convert all bases to uppercases\n");
		fprintf(stderr, "         -S        strip of white spaces in sequences\n");
		fprintf(stderr, "\n");
		free(kr);
		return 1;
	}
	if (line_len == 0) line_len = UINT_MAX;
	fp = optind < argc && strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	seq = kseq_init(fp);
	qual_thres += qual_shift;
	while (kseq_read(seq) >= 0) {
		++n_seqs;
		if (seq->seq.l < min_len) continue; // NB: length filter before taking random
		if (frac < 1. && kr_drand(kr) >= frac) continue;
		if (flag & 48) { // then choose odd/even reads only
			if ((flag&16) && (n_seqs&1) == 0) continue;
			if ((flag&32) && (n_seqs&1) == 1) continue;
		}
		if (flag & 512) { // option -S: squeeze out white spaces
			int k;
			if (seq->qual.l) {
				for (i = k = 0; i < seq->seq.l; ++i)
					if (!isspace(seq->seq.s[i]))
						seq->qual.s[k++] = seq->qual.s[i];
				seq->qual.l = k;
			}
			for (i = k = 0; i < seq->seq.l; ++i)
				if (!isspace(seq->seq.s[i]))
					seq->seq.s[k++] = seq->seq.s[i];
			seq->seq.l = k;
		}
		if (seq->qual.l && qual_thres > qual_shift) {
			if (mask_chr) {
				for (i = 0; i < seq->seq.l; ++i)
					if (seq->qual.s[i] < qual_thres || seq->qual.s[i] > max_q)
						seq->seq.s[i] = mask_chr;
			} else {
				for (i = 0; i < seq->seq.l; ++i)
					if (seq->qual.s[i] < qual_thres || seq->qual.s[i] > max_q)
						seq->seq.s[i] = tolower(seq->seq.s[i]);
			}
		}
		if (flag & 256) // option -U: convert to uppercases
			for (i = 0; i < seq->seq.l; ++i)
				seq->seq.s[i] = toupper(seq->seq.s[i]);
		if (flag & 1) seq->qual.l = 0; // option -a: fastq -> fasta
		if (flag & 2) seq->comment.l = 0; // option -C: drop fasta/q comments
		if (h) stk_mask(seq, h, flag&8, mask_chr); // masking
		if (flag & 4) { // option -r: reverse complement
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
		if (flag & 128) { // option -N: drop sequences containing ambiguous bases - Note: this is the last step!
			for (i = 0; i < seq->seq.l; ++i)
				if (seq_nt16to4_table[seq_nt16_table[(int)seq->seq.s[i]]] > 3) break;
			if (i < seq->seq.l) continue;
		}
		stk_printseq(seq, line_len);
	}
	kseq_destroy(seq);
	gzclose(fp);
	stk_reg_destroy(h);
	free(kr);
	return 0;
}

int stk_gc(int argc, char *argv[])
{
	int c, is_at = 0, min_l = 20;
	double frac = 0.6f, xdropoff = 10.0f, q;
	gzFile fp;
	kseq_t *seq;

	while ((c = getopt(argc, argv, "wx:f:l:")) >= 0) {
		if (c == 'x') xdropoff = atof(optarg);
		else if (c == 'w') is_at = 1;
		else if (c == 'f') frac = atof(optarg);
		else if (c == 'l') min_l = atoi(optarg);
	}
	if (optind + 1 > argc) {
		fprintf(stderr, "Usage: seqtk gc [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -w         identify high-AT regions\n");
		fprintf(stderr, "  -f FLOAT   min GC fraction (or AT fraction for -w) [%.2f]\n", frac);
		fprintf(stderr, "  -l INT     min region length to output [%d]\n", min_l);
		fprintf(stderr, "  -x FLOAT   X-dropoff [%.1f]\n", xdropoff);
		return 1;
	}
	q = (1.0f - frac) / frac;

	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	seq = kseq_init(fp);
	while (kseq_read(seq) >= 0) {
		int i, start = 0, max_i = 0, n_hits = 0, start_hits = 0, max_hits = 0;
		double sc = 0., max = 0.;
		for (i = 0; i < seq->seq.l; ++i) {
			int hit;
			c = seq_nt16_table[(int)seq->seq.s[i]];
			if (is_at) hit = (c == 1 || c == 8 || c == 9);
			else hit = (c == 2 || c == 4 || c == 6);
			n_hits += hit;
			if (hit) {
				if (sc == 0) start = i, start_hits = n_hits;
				sc += q;
				if (sc > max) max = sc, max_i = i, max_hits = n_hits;
			} else if (sc > 0) {
				sc += -1.0f;
				if (sc < 0 || max - sc > xdropoff) {
					if (max_i + 1 - start >= min_l)
						printf("%s\t%d\t%d\t%d\n", seq->name.s, start, max_i + 1, max_hits - start_hits + 1);
					sc = max = 0;
					i = max_i;
				}
			}
		}
		if (max > 0. && max_i + 1 - start >= min_l)
			printf("%s\t%d\t%d\t%d\n", seq->name.s, start, max_i + 1, max_hits - start_hits + 1);
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int stk_mergepe(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *seq[2];

	if (argc < 3) {
		fprintf(stderr, "Usage: seqtk mergepe <in1.fq> <in2.fq>\n");
		return 1;
	}
	fp1 = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	fp2 = strcmp(argv[2], "-")? gzopen(argv[2], "r") : gzdopen(fileno(stdin), "r");
	if (fp1 == 0 || fp2 == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	seq[0] = kseq_init(fp1);
	seq[1] = kseq_init(fp2);
	while (kseq_read(seq[0]) >= 0) {
		if (kseq_read(seq[1]) < 0) {
			fprintf(stderr, "[W::%s] the 2nd file has fewer records.\n", __func__);
			break;
		}
		stk_printseq(seq[0], 0);
		stk_printseq(seq[1], 0);
	}
	if (kseq_read(seq[1]) >= 0)
		fprintf(stderr, "[W::%s] the 1st file has fewer records.\n", __func__);
	kseq_destroy(seq[0]); gzclose(fp1);
	kseq_destroy(seq[1]); gzclose(fp2);
	return 0;
}

int stk_dropse(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq, last;

	if (argc == 1 && isatty(fileno(stdin))) {
		fprintf(stderr, "Usage: seqtk dropse <in.fq>\n");
		return 1;
	}
	fp = argc > 1 && strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	seq = kseq_init(fp);

	memset(&last, 0, sizeof(kseq_t));
	while (kseq_read(seq) >= 0) {
		if (last.name.l) {
			kstring_t *p = &last.name, *q = &seq->name;
			int is_diff;
			if (p->l == q->l) {
				int l = (p->l > 2 && p->s[p->l-2] == '/' && q->s[q->l-2] == '/' && isdigit(p->s[p->l-1]) && isdigit(q->s[q->l-1]))? p->l - 2 : p->l;
				is_diff = strncmp(p->s, q->s, l);
			} else is_diff = 1;
			if (!is_diff) {
				stk_printseq(&last, 0);
				stk_printseq(seq, 0);
				last.name.l = 0;
			} else cpy_kseq(&last, seq);
		} else cpy_kseq(&last, seq);
	}

	kseq_destroy(seq);
	gzclose(fp);
	// free last!
	return 0;
}

int stk_rename(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq, last;
	char *prefix = 0;
	uint64_t n = 1;

	if (argc == 1 && isatty(fileno(stdin))) {
		fprintf(stderr, "Usage: seqtk rename <in.fq> [prefix]\n");
		return 1;
	}
	fp = argc > 1 && strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	seq = kseq_init(fp);
	if (argc > 2) prefix = argv[2];

	memset(&last, 0, sizeof(kseq_t));
	while (kseq_read(seq) >= 0) {
		if (last.name.l) {
			kstring_t *p = &last.name, *q = &seq->name;
			int is_diff;
			if (p->l == q->l) {
				int l = (p->l > 2 && p->s[p->l-2] == '/' && q->s[q->l-2] == '/' && isdigit(p->s[p->l-1]) && isdigit(q->s[q->l-1]))? p->l - 2 : p->l;
				is_diff = strncmp(p->s, q->s, l);
			} else is_diff = 1;
			if (!is_diff) {
				stk_printseq_renamed(&last, 0, prefix, n);
				stk_printseq_renamed(seq,   0, prefix, n);
				last.name.l = 0;
				++n;
			} else {
				stk_printseq_renamed(&last, 0, prefix, n);
				++n;
				cpy_kseq(&last, seq);
			}
		} else cpy_kseq(&last, seq);
	}
	if (last.name.l) stk_printseq_renamed(&last, 0, prefix, n);

	kseq_destroy(seq);
	gzclose(fp);
	// free last!
	return 0;
}

int stk_kfreq(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *ks;
	int kmer, i, l, mask;
	char *nei;

	if (argc < 2) {
		fprintf(stderr, "Usage: seqtk kfreq <kmer> <in.fa>\n");
		return 1;
	}

	// get the k-mer
	l = strlen(argv[1]);
	for (i = kmer = 0; i < l; ++i) {
		int c = seq_nt6_table[(int)argv[1][i]];
		assert(c >= 1 && c <= 4);
		kmer = kmer << 2 | (c - 1);
	}
	mask = (1<<2*l) - 1;

	// get the neighbors
	nei = calloc(1, 1<<2*l);
	for (i = 0; i < l; ++i) {
		int j, x;
		x = kmer & ~(3 << 2*i);
		for (j = 0; j < 4; ++j) 
			nei[x|j<<2*i] = 1;
	}

	fp = argc == 2 || strcmp(argv[2], "-") == 0? gzdopen(fileno(stdin), "r") : gzopen(argv[2], "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		int k, x[2], cnt[2], cnt_nei[2], which;
		x[0] = x[1] = k = cnt[0] = cnt[1] = cnt_nei[0] = cnt_nei[1] = 0;
		for (i = 0; i < ks->seq.l; ++i) {
			int c = seq_nt6_table[(int)ks->seq.s[i]];
			if (c >= 1 && c <= 4) {
				x[0] = (x[0] << 2 | (c - 1)) & mask;
				x[1] = (x[1] >> 2 | (4 - c) << 2*(l-1));
				if (k < l) ++k;
				if (k == l) {
					if (x[0] == kmer) ++cnt[0];
					else if (x[1] == kmer) ++cnt[1];
					if (nei[x[0]]) ++cnt_nei[0];
					else if (nei[x[1]]) ++cnt_nei[1];
				}
			} else k = 0;
		}
		which = cnt_nei[0] > cnt_nei[1]? 0 : 1;
		printf("%s\t%ld\t%c\t%d\t%d\n", ks->name.s, ks->seq.l, "+-"[which], cnt_nei[which], cnt[which]);
	}
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}

/* fqchk */

typedef struct {
	int64_t q[94], b[5];
} posstat_t;

static void fqc_aux(posstat_t *p, int pos, int64_t allq[94], double perr[94], int qthres)
{
	int k;
	int64_t sum = 0, qsum = 0, sum_low = 0;
	double psum = 0;
	if (pos <= 0) printf("ALL");
	else printf("%d", pos);
	for (k = 0; k <= 4; ++k) sum += p->b[k];
	printf("\t%lld", (long long)sum);
	for (k = 0; k <= 4; ++k)
		printf("\t%.1f", 100. * p->b[k] / sum);
	for (k = 0; k <= 93; ++k) {
		qsum += p->q[k] * k, psum += p->q[k] * perr[k];
		if (k < qthres) sum_low += p->q[k];
	}
	printf("\t%.1f\t%.1f", (double)qsum/sum, -4.343*log((psum+1e-6)/(sum+1e-6)));
	if (qthres <= 0) {
		for (k = 0; k <= 93; ++k)
			if (allq[k] > 0) printf("\t%.2f", 100. * p->q[k] / sum);
	} else printf("\t%.1f\t%.1f", 100. * sum_low / sum, 100. * (sum - sum_low) / sum);
	putchar('\n');
}

int stk_fqchk(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int i, c, k, max_len = 0, min_len = 0x7fffffff, max_alloc = 0, offset = 33, n_diffQ = 0, qthres = 20;
	int64_t tot_len = 0, n = 0;
	double perr[94];
	posstat_t all, *pos = 0;

	while ((c = getopt(argc, argv, "q:")) >= 0)
		if (c == 'q') qthres = atoi(optarg);

	if (optind == argc) {
		fprintf(stderr, "Usage: seqtk fqchk [-q %d] <in.fq>\n", qthres);
		fprintf(stderr, "Note: use -q0 to get the distribution of all quality values\n");
		return 1;
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
		return 1;
	}
	seq = kseq_init(fp);
	for (k = 0; k <= 93; ++k)
		perr[k] = pow(10., -.1 * k);
	perr[0] = perr[1] = perr[2] = perr[3] = .5;
	while (kseq_read(seq) >= 0) {
		if (seq->qual.l == 0) continue;
		++n;
		tot_len += seq->seq.l;
		min_len = min_len < seq->seq.l? min_len : seq->seq.l;
		max_len = max_len > seq->seq.l? max_len : seq->seq.l;
		if (max_len > max_alloc) {
			int old_max = max_alloc;
			max_alloc = max_len;
			kroundup32(max_alloc);
			pos = realloc(pos, max_alloc * sizeof(posstat_t));
			memset(&pos[old_max], 0, (max_alloc - old_max) * sizeof(posstat_t));
		}
		for (i = 0; i < seq->qual.l; ++i) {
			int q = seq->qual.s[i] - offset;
			int b = seq_nt6_table[(int)seq->seq.s[i]];
			b = b? b - 1 : 4;
			q = q < 93? q : 93;
			++pos[i].q[q];
			++pos[i].b[b];
		}
	}
	kseq_destroy(seq);
	gzclose(fp);

	memset(&all, 0, sizeof(posstat_t));
	for (i = 0; i < max_len; ++i) {
		for (k = 0; k <= 93; ++k)
			all.q[k] += pos[i].q[k];
		for (k = 0; k <= 4; ++k)
			all.b[k] += pos[i].b[k];
	}
	for (k = n_diffQ = 0; k <= 93; ++k)
		if (all.q[k]) ++n_diffQ;
	printf("min_len: %d; max_len: %d; avg_len: %.2f; %d distinct quality values\n", min_len, max_len, (double)tot_len/n, n_diffQ);
	printf("POS\t#bases\t%%A\t%%C\t%%G\t%%T\t%%N\tavgQ\terrQ");
	if (qthres <= 0) {
		for (k = 0; k <= 93; ++k)
			if (all.q[k] > 0) printf("\t%%Q%d", k);
	} else printf("\t%%low\t%%high");
	putchar('\n');
	fqc_aux(&all, 0, all.q, perr, qthres);
	for (i = 0; i < max_len; ++i)
		fqc_aux(&pos[i], i + 1, all.q, perr, qthres);
	free(pos);
	return 0;
}

/* main function */
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   seqtk <command> <arguments>\n");
	fprintf(stderr, "Version: 1.2-r101-dirty\n\n");
	fprintf(stderr, "Command: seq       common transformation of FASTA/Q\n");
	fprintf(stderr, "         comp      get the nucleotide composition of FASTA/Q\n");
	fprintf(stderr, "         sample    subsample sequences\n");
	fprintf(stderr, "         subseq    extract subsequences from FASTA/Q\n");
	fprintf(stderr, "         fqchk     fastq QC (base/quality summary)\n");
	fprintf(stderr, "         mergepe   interleave two PE FASTA/Q files\n");
	fprintf(stderr, "         trimfq    trim FASTQ using the Phred algorithm\n\n");
	fprintf(stderr, "         hety      regional heterozygosity\n");
	fprintf(stderr, "         gc        identify high- or low-GC regions\n");
	fprintf(stderr, "         mutfa     point mutate FASTA at specified positions\n");
	fprintf(stderr, "         mergefa   merge two FASTA/Q files\n");
	fprintf(stderr, "         famask    apply a X-coded FASTA to a source FASTA\n");
	fprintf(stderr, "         dropse    drop unpaired from interleaved PE FASTA/Q\n");
	fprintf(stderr, "         rename    rename sequence names\n");
	fprintf(stderr, "         randbase  choose a random base from hets\n");
	fprintf(stderr, "         cutN      cut sequence at long N\n");
	fprintf(stderr, "         listhet   extract the position of each het\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();
	if (strcmp(argv[1], "comp") == 0) stk_comp(argc-1, argv+1);
	else if (strcmp(argv[1], "fqchk") == 0) stk_fqchk(argc-1, argv+1);
	else if (strcmp(argv[1], "hety") == 0) stk_hety(argc-1, argv+1);
	else if (strcmp(argv[1], "gc") == 0) stk_gc(argc-1, argv+1);
	else if (strcmp(argv[1], "subseq") == 0) stk_subseq(argc-1, argv+1);
	else if (strcmp(argv[1], "mutfa") == 0) stk_mutfa(argc-1, argv+1);
	else if (strcmp(argv[1], "mergefa") == 0) stk_mergefa(argc-1, argv+1);
	else if (strcmp(argv[1], "mergepe") == 0) stk_mergepe(argc-1, argv+1);
	else if (strcmp(argv[1], "dropse") == 0) stk_dropse(argc-1, argv+1);
	else if (strcmp(argv[1], "randbase") == 0) stk_randbase(argc-1, argv+1);
	else if (strcmp(argv[1], "cutN") == 0) stk_cutN(argc-1, argv+1);
	else if (strcmp(argv[1], "listhet") == 0) stk_listhet(argc-1, argv+1);
	else if (strcmp(argv[1], "famask") == 0) stk_famask(argc-1, argv+1);
	else if (strcmp(argv[1], "trimfq") == 0) stk_trimfq(argc-1, argv+1);
	else if (strcmp(argv[1], "hrun") == 0) stk_hrun(argc-1, argv+1);
	else if (strcmp(argv[1], "sample") == 0) stk_sample(argc-1, argv+1);
	else if (strcmp(argv[1], "seq") == 0) stk_seq(argc-1, argv+1);
	else if (strcmp(argv[1], "kfreq") == 0) stk_kfreq(argc-1, argv+1);
	else if (strcmp(argv[1], "rename") == 0) stk_rename(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
