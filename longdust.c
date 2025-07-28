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

static double ld_f_large(double lambda) // with Sterling's approximation
{
	double x = 0.5 * log(2.0 * M_PI * M_E * lambda) - 1.0 / 12.0 / lambda * (1.0 + 0.5 / lambda + 19.0 / 30.0 / lambda / lambda);
	x += lambda * (log(lambda) - 1.0);
	return x;
}

static double *ld_cal_f(void *km, int32_t k, int32_t max_l)
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

static unsigned char seq_comp_tab[] = { // not really necessary
	  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
	 16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
	 32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
	 48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

/****************************
 * Internal data structures *
 ****************************/

KDQ_INIT(uint32_t)

struct ld_data_s {
	void *km; // memory
	const ld_opt_t *opt;
	double *f, *c;
	kdq_t(uint32_t) *q;
	int32_t *ht;
	// output
	int64_t n_intv, m_intv;
	ld_intv_t *intv;
};

void ld_opt_init(ld_opt_t *opt)
{
	opt->kmer = 6;
	opt->ws = 1024;
	opt->thres = 0.6;
	opt->xdrop_len1 = 200;
	opt->xdrop_len2 = 50;
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

static int32_t ld_dust_back(ld_data_t *ld, int64_t pos, double ws_sum) // pos only for debugging
{
	const ld_opt_t *opt = ld->opt;
	double xdrop1 = opt->thres * opt->xdrop_len1;
	double xdrop2 = opt->thres * opt->xdrop_len2;
	int32_t i, l, max_i;
	double s, sl, max_sf = 0.0, max_sb = 0.0;

	// backward
	memset(ld->ht, 0, sizeof(int32_t) * (1U<<2*opt->kmer));
	for (i = kdq_size(ld->q) - 1, l = 1, s = sl = 0.0, max_i = -1; i >= 0; --i, ++l) {
		uint32_t x = kdq_at(ld->q, i);
		s += (x&1? 0 : ld->c[++ld->ht[x>>1]]) - opt->thres;
		sl = s - ld->f[l];
		if (sl > max_sb) {
			max_sb = sl, max_i = i;
		} else if (max_sb == 0.0) { // haven't gone beyond the baseline before
			if (max_sb - sl > xdrop1) break;
		} else { // max_sb > 0.0 in this case
			if (max_sb - sl > xdrop2) break;
		}
		if (ws_sum - ld->f[l] - l * opt->thres < max_sb) break;
	}
	if (max_i < 0) return -1;

	// forward
	memset(ld->ht, 0, sizeof(int32_t) * (1U<<2*opt->kmer));
	for (i = max_i, l = 1, s = sl = 0.0; i < kdq_size(ld->q); ++i, ++l) {
		uint32_t x = kdq_at(ld->q, i);
		s += (x&1? 0 : ld->c[++ld->ht[x>>1]]) - opt->thres;
		sl = s - ld->f[l];
		if (sl >= max_sf) max_sf = sl;
	}
	//fprintf(stderr, "[%ld,%ld]:%c\t(%f,%f,%f)\n", (long)(pos-(kdq_size(ld->q)-1-max_i)-(opt->kmer-1)), (long)pos, "ACGT"[kdq_at(ld->q, kdq_size(ld->q)-1)>>1&3], max_sb, max_sf, sl);
	return sl >= max_sf - 1e-6? max_i : -1;
}

static int32_t ld_extend(ld_data_t *ld)
{
	uint32_t x = kdq_at(ld->q, kdq_size(ld->q) - 1);
	int32_t l = kdq_size(ld->q) - 1;
	double diff;
	if (x&1) return -1;
	diff = ld->c[ld->ht[x>>1] + 1] - (ld->f[l + 1] - ld->f[l]);
	if (diff < ld->opt->thres) return -1; // if this doesn't increase score, don't extend
	return 0;
}

static void ld_save_intv(ld_data_t *ld, int64_t st, int64_t en)
{
	int64_t k;
	for (k = ld->n_intv - 1; k >= 0; --k) // intv[] is sorted by end positions
		if (st > ld->intv[k].en)
			break;
	++k;
	if (k < ld->n_intv) { // overlapping with one (or multiple) saved interval; overwrite the leftmost saved one
		if (ld->intv[k].st > st) ld->intv[k].st = st;
		if (ld->intv[k].en < en) ld->intv[k].en = en;
		ld->n_intv = k + 1;
	} else { // create a new entry
		Kgrow(ld->km, ld_intv_t, ld->intv, ld->n_intv, ld->m_intv);
		ld->intv[ld->n_intv].st = st;
		ld->intv[ld->n_intv++].en = en;
	}
}

void ld_dust1(ld_data_t *ld, int64_t len, const uint8_t *seq)
{
	const ld_opt_t *opt = ld->opt;
	uint32_t x, mask = (1U<<2*opt->kmer) - 1;
	int64_t i, l, st = -1, en = -1;
	int64_t last_i = -1, last_q = -1; // last_i: i of last successful LCR; last_q: pos in queue when last_i is set
	int32_t *ht;
	double ht_sum = 0.0;

	ld->n_intv = 0;
	ld->q->front = ld->q->count = 0;
	ht = Kcalloc(ld->km, int32_t, mask + 1);
	for (i = 0, x = 0, l = 0; i <= len; ++i) {
		int32_t ambi, j, b = i < len? seq_nt4_table[seq[i]] : 4;
		if (b < 4) {
			x = (x << 2 | b) & mask;
			++l;
			ambi = (l < opt->kmer);
		} else {
			l = 0;
			ambi = 1;
		}
		if (kdq_size(ld->q) >= opt->ws) { // remove from the queue
			uint32_t *p;
			p = kdq_shift(uint32_t, ld->q);
			if ((*p&1) == 0) ht_sum -= ld->c[ht[*p>>1]--];
			if (last_q == 0) --ld->ht[*p>>1];
			else --last_q; // this needs to be updated as the queue is shifted
		}
		kdq_push(uint32_t, ld->q, x<<1|ambi);
		if (ambi) continue;

		j = -1;
		ht_sum += ld->c[++ht[x]];
		if (ht[x] > 1) { // no need to call ld_dust_back() if x is a singleton in the window
			if (last_i == i - 1 && last_q == 0) j = ld_extend(ld); // test and potentially extend the base at i
			if (j < 0) j = ld_dust_back(ld, i, ht_sum); // do full dust_back
		}
		if (j >= 0) { // LCR found
			int64_t st2 = i - (kdq_size(ld->q) - 1 - j) - (opt->kmer - 1); // the start of LCR
			if (st2 < en) { // overlapping with the active LCR interval
				if (st < 0 || st2 < st) st = st2;
			} else { // not overlapping; save it to intv[]
				if (st >= 0) ld_save_intv(ld, st, en);
				st = st2;
			}
			en = i + 1;
			last_i = i, last_q = j;
		}
	}
	if (st >= 0) ld_save_intv(ld, st, en);
	kfree(ld->km, ht);
}

void ld_dust2(ld_data_t *ld, int64_t len, const uint8_t *seq)
{
	int64_t i, j[2], st, en, n[2];
	uint8_t *rev;
	ld_intv_t *intv[2];

	// forward
	ld_dust1(ld, len, seq);
	n[0] = ld->n_intv;
	intv[0] = Kmalloc(ld->km, ld_intv_t, n[0]);
	memcpy(intv[0], ld->intv, n[0] * sizeof(ld_intv_t));

	// reverse
	rev = Kmalloc(ld->km, uint8_t, len);
	for (i = 0; i < len; ++i)
		rev[len - i - 1] = seq[i] < 128? seq_comp_tab[seq[i]] : seq[i];
	ld_dust1(ld, len, rev);
	kfree(ld->km, rev);
	n[1] = ld->n_intv;
	intv[1] = Kmalloc(ld->km, ld_intv_t, n[1]);
	for (i = 0; i < ld->n_intv; ++i) {
		intv[1][ld->n_intv - 1 - i].st = len - ld->intv[i].en;
		intv[1][ld->n_intv - 1 - i].en = len - ld->intv[i].st;
	}

	// merge
	ld->n_intv = 0;
	st = en = 0;
	j[0] = j[1] = 0;
	while (j[0] < n[0] || j[1] < n[1]) {
		int32_t w = j[0] >= n[0]? 1 : j[1] >= n[1]? 0 : intv[0][j[0]].st < intv[1][j[1]].st? 0 : 1;
		const ld_intv_t *p = &intv[w][j[w]++];
		if (p->st <= en) {
			en = en > p->en? en : p->en;
		} else {
			if (en > st) {
				Kgrow(ld->km, ld_intv_t, ld->intv, ld->n_intv, ld->m_intv);
				ld->intv[ld->n_intv].st = st;
				ld->intv[ld->n_intv++].en = en;
			}
			st = p->st, en = p->en;
		}
	}
	if (en > st) {
		Kgrow(ld->km, ld_intv_t, ld->intv, ld->n_intv, ld->m_intv);
		ld->intv[ld->n_intv].st = st;
		ld->intv[ld->n_intv++].en = en;
	}
	kfree(ld->km, intv[1]);
	kfree(ld->km, intv[0]);
}

const ld_intv_t *ld_get_intv(const ld_data_t *ld, int64_t *n)
{
	*n = ld->n_intv;
	return ld->intv;
}
