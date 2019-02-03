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

#ifndef __SEQTK_H__
#define __SEQTK_H__

#include "khash.h"
#include "constants.h"
#include <inttypes.h>
#include "kseq.h"

typedef struct {
	int n, m;
	uint64_t *a;
} reglist_t;

KHASH_MAP_INIT_STR(reg, reglist_t)
KHASH_SET_INIT_INT64(64)

typedef kh_reg_t reghash_t;

typedef uint64_t krint64_t;

struct _krand_t;
typedef struct _krand_t krand_t;

struct _krand_t {
    int mti;
    krint64_t mt[KR_NN];
};

typedef struct {
    int64_t q[94], b[5];
} posstat_t;

extern reghash_t *stk_reg_read(const char *fn);
extern void       stk_reg_destroy(reghash_t *h);
extern krand_t   *kr_srand(krint64_t seed);
extern krint64_t  kr_rand(krand_t *kr);
extern int        stk_trimfq(int argc, char *argv[]);
extern int        stk_comp(int argc, char *argv[]);
extern int        stk_randbase(int argc, char *argv[]);
extern int        stk_hety(int argc, char *argv[]);
extern int        stk_subseq(int argc, char *argv[]);
extern int        stk_mergefa(int argc, char *argv[]);
extern int        stk_famask(int argc, char *argv[]);
extern int        stk_mutfa(int argc, char *argv[]);
extern int        stk_listhet(int argc, char *argv[]);
extern int        stk_cutN(int argc, char *argv[]);
extern int        stk_hrun(int argc, char *argv[]);
extern int        stk_sample(int argc, char *argv[]);
extern int        stk_seq(int argc, char *argv[]);
extern int        stk_gc(int argc, char *argv[]);
extern int        stk_mergepe(int argc, char *argv[]);
extern int        stk_dropse(int argc, char *argv[]);
extern int        stk_rename(int argc, char *argv[]);
extern int        stk_kfreq(int argc, char *argv[]);
extern int        stk_fqchk(int argc, char *argv[]);


#endif
