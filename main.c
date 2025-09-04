#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "longdust.h"
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c, for_only = 0;
	gzFile fp;
	kseq_t *ks;
	ld_opt_t opt;
	ld_data_t *ld;

	ld_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "k:w:e:t:fva", 0)) >= 0) {
		if (c == 'k') opt.kmer = atoi(o.arg);
		else if (c == 'w') opt.ws = atoi(o.arg);
		else if (c == 't') opt.thres = atof(o.arg);
		else if (c == 'e') opt.xdrop_len = atoi(o.arg);
		else if (c == 'a') opt.approx = 1;
		else if (c == 'f') for_only = 1;
		else if (c == 'v') {
			puts(LD_VERSION);
			return 0;
		}
	}
	if (argc - o.ind == 0) {
		fprintf(stderr, "Usage: longdust [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT      k-mer length [%d]\n", opt.kmer);
		fprintf(stderr, "  -w INT      window size [%d]\n", opt.ws);
		fprintf(stderr, "  -t FLOAT    score threshold [%g]\n", opt.thres);
		fprintf(stderr, "  -e INT      extension X-drop length (0 to disable) [%d]\n", opt.xdrop_len);
		fprintf(stderr, "  -f          forward strand only\n");
		fprintf(stderr, "  -a          guaranteed O(Lw) algorithm but with more approximation\n");
		fprintf(stderr, "  -v          version number\n");
		fprintf(stderr, "Notes:\n");
		fprintf(stderr, "  * Recommend w < 4^k for performance, especially given large w\n");
		fprintf(stderr, "  * Use \"-k6 -w1000 -t.55\" for more relaxed but shorter regions\n");
		return 1;
	}

	fp = strcmp(argv[o.ind], "-")? gzopen(argv[o.ind], "r") : gzdopen(0, "r");
	ks = kseq_init(fp);
	ld = ld_data_init(0, &opt);
	while (kseq_read(ks) >= 0) {
		int64_t i, n;
		const ld_intv_t *intv;
		if (for_only) ld_dust1(ld, ks->seq.l, (uint8_t*)ks->seq.s);
		else ld_dust2(ld, ks->seq.l, (uint8_t*)ks->seq.s);
		intv = ld_get_intv(ld, &n);
		for (i = 0; i < n; ++i)
			printf("%s\t%ld\t%ld\n", ks->name.s, (long)intv[i].st, (long)intv[i].en);
	}
	ld_data_destroy(ld);
	kseq_destroy(ks);
	return 0;
}
