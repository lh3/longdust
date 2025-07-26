#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "kalloc.h"
#include "kdq.h"
#include "longdust.h"

/***************
 * Compute f() *
 ***************/

double ld_f_large(double lambda) // with Sterling's approximation
{
	double x = 0.5 * log(2.0 * M_PI * M_E * lambda) - 1.0 / 12.0 / lambda * (1.0 + 0.5 / lambda + 19.0 / 30.0 / lambda / lambda);
	x += lambda * (log(lambda) - 1.0);
	return x;
}

double *ld_cal_f(void *km, int32_t k, int32_t max_l)
{
	static const double eps = 1e-9;
	static const int32_t max_n = 10000;
	double *f;
	int32_t l;
	assert(k < 16);
	f = Kcalloc(km, double, max_l + 1);
	for (l = 1; l <= max_l; ++l) {
		double lambda = (double)l / (1U<<2*k);
		double x = 0.0, sn = 0.0, y = lambda;
		int32_t n;
		if (lambda < 30.0) {
			for (n = 2; n <= max_n; ++n) {
				double z;
				sn += log(n);
				y *= lambda / n;
				z = y * sn;
				if (z < x * eps) break;
				x += z;
			}
			f[l] = x * exp(-lambda);
			//printf("%d\t%f\t%f\t%f\n", l, lambda, f[l], ld_f_large(lambda));
		} else {
			f[l] = ld_f_large(lambda);
		}
	}
	return f;
}

/*********
 * Table *
 *********/

static unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
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

/****************************
 * Internal data structures *
 ****************************/

KDQ_INIT(uint32_t)

struct ld_data_s {
	const ld_opt_t *opt;
	void *km;
	kdq_t(uint32_t) *q;
	int32_t *ht;
	double *f, *c;
	int64_t n_intv, m_intv;
	ld_intv_t *intv;
};

void ld_opt_init(ld_opt_t *opt)
{
	opt->kmer = 6;
	opt->ws = 1024;
	opt->thres = 0.5;
	opt->xdrop = 100.0;
}

ld_data_t *ld_data_init(void *km, const ld_opt_t *opt)
{
	int32_t i;
	ld_data_t *ld;
	ld = Kcalloc(km, ld_data_t, 1);
	ld->opt = opt;
	ld->ht = Kcalloc(km, int32_t, 1U<<2*opt->kmer);
	ld->f = ld_cal_f(km, opt->kmer, opt->ws);
	ld->c = Kcalloc(km, double, opt->ws + 1);
	for (i = 2; i <= opt->ws; ++i)
		ld->c[i] = log(i);
	ld->q = kdq_init(uint32_t, km);
	return ld;
}

void ld_data_destroy(ld_data_t *ld)
{
	void *km = ld->km;
	kfree(km, ld->intv);
	kdq_destroy(uint32_t, ld->q);
	kfree(km, ld->c); kfree(km, ld->f);
	kfree(km, ld->ht);
	kfree(km, ld);
}

int32_t ld_dust_pos(ld_data_t *ld)
{
	const ld_opt_t *opt = ld->opt;
	int32_t i, max_i, cnt;
	double s, max_sf, max_sb;

	// backward
	memset(ld->ht, 0, sizeof(int32_t) * (1U<<2*opt->kmer));
	for (i = kdq_size(ld->q) - 1, s = 0.0, max_sb = 0.0, max_i = -1; i >= 0; --i) {
		uint32_t x;
		x = kdq_at(ld->q, i);
		if ((x & 1) == 0) {
			cnt = ++ld->ht[x>>1];
			//if (kdq_size(ld->q) == 512) fprintf(stderr, "X2\t%d\t%x\t%d\t%f\n", i, x, cnt, s);
			s += ld->c[cnt] - opt->thres;
		} else {
			s -= opt->thres;
		}
		if (s > max_sb)
			max_sb = s, max_i = i;
		else if (max_sb - s > opt->xdrop)
			break;
	}
	if (max_i < 0) return -1;

	// forward
	memset(ld->ht, 0, sizeof(int32_t) * (1U<<2*opt->kmer));
	for (i = max_i, s = 0.0, max_sf = 0.0; i < kdq_size(ld->q); ++i) {
		uint32_t x;
		x = kdq_at(ld->q, i);
		if ((x & 1) == 0) {
			cnt = ++ld->ht[x>>1]; // this is the new count
			s += ld->c[cnt] - opt->thres;
		} else {
			s -= opt->thres;
		}
		if (s >= max_sf) max_sf = s;
	}
	//if (kdq_size(ld->q) == 512) fprintf(stderr, "X1\t%d\t%f\t%f\t%f\n", max_i, max_sb, max_sf, s);
	return s >= max_sf - 1e-6? max_i : -1;
}

void ld_dust(ld_data_t *ld, int64_t len, const uint8_t *seq)
{
	const ld_opt_t *opt = ld->opt;
	uint32_t x, mask = (1U<<2*opt->kmer) - 1;
	int64_t i, l, st = -1, en = -1;
	ld->n_intv = 0;
	for (i = 0, x = 0, l = 0; i <= len; ++i) {
		int32_t j, ambi, b = i < len? seq_nt4_table[seq[i]] : 4;
		if (i%1000000 == 0) fprintf(stderr, "%ld\t%ld\n", (long)i, (long)kdq_size(ld->q));
		if (b < 4) {
			x = (x << 2 | b) & mask;
			++l;
			ambi = (l < opt->kmer);
		} else {
			l = 0;
			ambi = 1;
		}
		if (kdq_size(ld->q) >= opt->ws)
			kdq_shift(uint32_t, ld->q);
		kdq_push(uint32_t, ld->q, x<<1|ambi);
		if (!ambi) {
			j = ld_dust_pos(ld);
			if (j > 0) {
				int64_t st2 = i - (kdq_size(ld->q) - 1 - j) - (opt->kmer - 1);
				//if (kdq_size(ld->q) == 512) fprintf(stderr, "[%ld,%ld)\t%ld\t[%ld,%ld)\n", (long)st2, (long)i, (long)j, (long)st, (long)en);
				if (st2 < en) {
					if (st < 0 || st2 < st) st = st2;
				} else {
					if (st >= 0 && en - st >= opt->kmer) {
						Kgrow(ld->km, ld_intv_t, ld->intv, ld->n_intv, ld->m_intv);
						ld->intv[ld->n_intv].st = st;
						ld->intv[ld->n_intv++].en = en;
					}
					st = st2;
				}
				en = i + 1;
			}
		}
	}
}

/*****************
 * main function *
 *****************/

#ifndef LD_LIB_ONLY
#include <zlib.h>
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int32_t c;
	gzFile fp;
	kseq_t *ks;
	ld_opt_t opt;
	ld_data_t *ld;

	ld_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "k:w:", 0)) >= 0) {
		if (c == 'k') opt.kmer = atoi(o.arg);
		else if (c == 'w') opt.ws = atoi(o.arg);
	}
	if (argc - o.ind == 0) {
		fprintf(stderr, "Usage: longdust [options] <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT      k-mer length [%d]\n", opt.kmer);
		fprintf(stderr, "  -w INT      window size [%d]\n", opt.ws);
		return 1;
	}

	fp = strcmp(argv[o.ind], "-")? gzopen(argv[o.ind], "r") : gzdopen(0, "r");
	ks = kseq_init(fp);
	ld = ld_data_init(0, &opt);
	while (kseq_read(ks) >= 0) {
		int64_t i;
		ld_dust(ld, ks->seq.l, (uint8_t*)ks->seq.s);
		for (i = 0; i < ld->n_intv; ++i) {
			printf("%s\t%ld\t%ld\n", ks->name.s, (long)ld->intv[i].st, (long)ld->intv[i].en);
		}
	}
	ld_data_destroy(ld);
	kseq_destroy(ks);
	return 0;
}
#endif
