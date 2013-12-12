#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include "kvec.h"
#include "kstring.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 65536)

typedef kvec_t(uint32_t) vec32_t;
typedef kvec_t(uint64_t) vec64_t;

int main_cut(int argc, char *argv[])
{
	vec64_t cols = {0,0,0}, buf = {0,0,0};
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0}, out = {0,0,0};
	int dret, sep = '\t', c, reorder = 0;

	while ((c = getopt(argc, argv, "rd:f:")) >= 0) {
		if (c == 'r') reorder = 1;
		else if (c == 'd') {
			if (strcmp(optarg, "isspace") == 0) sep = 256;
			else if (strlen(optarg) == 1) sep = optarg[0];
			else {
				fprintf(stderr, "[E::%s] invalid delimitor\n", __func__);
				return 1;
			}
		} else if (c == 'f') {
			int32_t beg, end, x;
			char *p = optarg;
			while ((x = strtol(p, &p, 10)) > 0) { // parse the field string
				beg = end = x;
				if (*p == '-')
					end = (x = strtol(p + 1, &p, 10)) >= beg? x : INT32_MAX;
				kv_push(uint64_t, cols, (uint64_t)(beg-1)<<32 | end);
				if (*p) ++p; // skip ','
				else break;
			}
			if (*p) {
				fprintf(stderr, "[E::%s] malformated field string. Abort!\n", __func__);
				return 1;
			}
		}
	}
	
	if (argc == optind && isatty(fileno(stdin))) {
		fprintf(stderr, "\nUsage: tabtk cut [options] [file.txt]\n\n");
		fprintf(stderr, "Options: -d CHAR     delimitor, a single CHAR or 'isspace' for both SPACE and TAB [TAB]\n");
		fprintf(stderr, "         -f STR      fields to cut; format identical to Unix cut [null]\n");
		fprintf(stderr, "         -r          reorder fields\n\n");
		return 1;
	}

	if (cols.n == 0) {
		fprintf(stderr, "[E::%s] no list of fields is specified.\n", __func__);
		return 2;
	}

	if (!reorder) {
		uint64_t *i;
		int j, k;
		for (i = cols.a + 1; i < cols.a + cols.n; ++i) // insertion sort
			if (*i < *(i - 1)) {
				uint64_t *j, tmp = *i;
				for (j = i; j > cols.a && tmp < *(j-1); --j) *j = *(j - 1);
				*j = tmp;
			}
		for (j = k = 1; j < cols.n; ++j) { // merge overlapping col regions
			if (cols.a[j]>>32 <= (uint32_t)cols.a[k-1])
				cols.a[k-1] = cols.a[k-1]>>32<<32 | (uint32_t)cols.a[j];
			else cols.a[k++] = cols.a[j];
		}
		cols.n = k;
	}

	fp = (optind == argc && !isatty(fileno(stdin))) || strcmp(argv[optind], "-") == 0? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	ks = ks_init(fp);
	while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
		int b, i;
		buf.n = 0; out.l = 0;
		for (i = b = 0; i <= str.l; ++i) { // mark columns
			if (str.s[i] == sep || i == str.l || (sep == 256 && isspace(str.s[i]))) {
				kv_push(uint64_t, buf, (uint64_t)b<<32 | i);
				b = i + 1;
			}
		}
		for (i = 0; i < cols.n; ++i) { // print columns
			int32_t j, beg = cols.a[i]>>32, end = (int32_t)cols.a[i];
			for (j = beg; j < end && j < buf.n; ++j) {
				uint64_t x = buf.a[j];
				if (out.l) kputc('\t', &out);
				kputsn(&str.s[x>>32], (uint32_t)x - (x>>32), &out);
			}
		}
		puts(out.s? out.s : "");
	}
	ks_destroy(ks);
	gzclose(fp);

	free(str.s); free(out.s); free(cols.a); free(buf.a);
	return 0;
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   tabtk <command> [arguments]\n");
	fprintf(stderr, "Version: 0.0\n\n");
	fprintf(stderr, "Command: cut       Unix cut with optional column reordering\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();
	if (strcmp(argv[1], "cut") == 0) main_cut(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized commad '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
