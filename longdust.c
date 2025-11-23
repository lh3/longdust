#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "kalloc.h"
#include "kdq.h"
#include "longdust.h"

#define LD_MAX_K 15

/***************
 * Compute f() *
 ***************/

#define LD_E  2.71828182845904523536028747135266250
#define LD_PI 3.14159265358979323846264338327950288

typedef struct {
	int32_t k, dummy;
	double gc;
	double dr[LD_MAX_K + 1];
	int32_t n_dr[LD_MAX_K + 1];
} ld_raux_t;

static ld_raux_t *ld_raux_init(void *km, int32_t k, double gc)
{
	ld_raux_t *r;
	uint32_t n_kmer = 1U << 2*k, x;
	int32_t i;
	assert(k < LD_MAX_K);
	r = Kcalloc(km, ld_raux_t, 1);
	for (i = 0; i <= k; ++i)
		r->dr[i] = pow(gc / 0.5, i) * pow((1.0 - gc) / 0.5, k - i);
	for (x = 0; x < n_kmer; ++x) { // traverse each k-mer
		int32_t n_gc = 0;
		for (i = 0; i < k; ++i) // count GC on the k-mer
			if ((x>>2*i&3) == 1 || (x>>2*i&3) == 2)
				++n_gc;
		r->n_dr[n_gc]++;
	}
	return r;
}

static double ld_f_large(double lambda) // with Sterling's approximation
{
	double x = 0.5 * log(2.0 * LD_PI * LD_E * lambda) - 1.0 / 12.0 / lambda * (1.0 + 0.5 / lambda + 19.0 / 30.0 / lambda / lambda);
	return x + lambda * (log(lambda) - 1.0);
}

static double *ld_cal_f2(void *km, int32_t k, int32_t max_l, int32_t nn_dr, const int32_t *n_dr, const double *dr)
{
	static const double eps = 1e-9;
	static const int32_t max_n = 10000;
	double *f;
	int32_t i, l;
	uint32_t n_kmer = 1U << 2*k;

	assert(k < LD_MAX_K);
	f = Kcalloc(km, double, max_l + 1);
	for (l = 1; l <= max_l; ++l) {
		for (i = 0; i < nn_dr; ++i) {
			double lambda = (double)l / n_kmer * dr[i];
			double x = 0.0, sn = 0.0, y = lambda, fli;
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
				fli = x * exp(-lambda);
				//printf("%d\t%f\t%f\t%f\n", l, lambda, f[l], ld_f_large(lambda));
			} else {
				fli = ld_f_large(lambda);
			}
			f[l] += fli * n_dr[i];
		}
		//printf("%d\t%lf\n", l, f[l]);
	}
	return f;
}

static double *ld_cal_f(void *km, int32_t k, int32_t max_l)
{
	int32_t n_dr = 1U << 2*k;
	double dr = 1.0;
	return ld_cal_f2(km, k, max_l, 1, &n_dr, &dr);
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

typedef struct {
	int32_t pos;
	double max;
} ld_forpos_t;

struct ld_data_s {
	void *km; // memory
	const ld_opt_t *opt;
	double *f, *c;
	kdq_t(uint32_t) *q;
	uint16_t *ht;
	int32_t max_test, n_for_pos;
	ld_forpos_t *for_pos;
	// output
	int64_t n_intv, m_intv;
	ld_intv_t *intv;
};

void ld_opt_init(ld_opt_t *opt)
{
	opt->kmer = 7;
	opt->ws = 5000;
	opt->thres = 0.6;
	opt->xdrop_len = 50;
	opt->min_start_cnt = 3;
	opt->approx = 0;
	opt->gc = -1.0; // GC correction disabled by default
}

ld_data_t *ld_data_init(void *km, const ld_opt_t *opt)
{
	int32_t i;
	double s, sl;
	ld_data_t *ld;
	assert(opt->ws > 0 && opt->ws < 0xffff);
	ld = Kcalloc(km, ld_data_t, 1);
	ld->opt = opt;
	ld->ht = Kcalloc(km, uint16_t, 1U<<2*opt->kmer);
	if (opt->gc > 0.0 && opt->gc < 1.0) {
		ld_raux_t *r;
		r = ld_raux_init(km, opt->kmer, opt->gc);
		ld->f = ld_cal_f2(km, opt->kmer, opt->ws + 1, opt->kmer + 1, r->n_dr, r->dr);
		kfree(km, r);
	} else {
		ld->f = ld_cal_f(km, opt->kmer, opt->ws + 1);
	}
	ld->c = Kcalloc(km, double, opt->ws + 1);
	for (i = 2; i <= opt->ws; ++i)
		ld->c[i] = log(i);
	ld->q = kdq_init(uint32_t, km);
	ld->for_pos = Kcalloc(km, ld_forpos_t, opt->ws + 1);
	// calculate min_test, the max step used in ld_if_backward()
	for (i = 1, s = 0.0; i < opt->ws; ++i) { // i <- minimum j such that \sum_{c=1}^j {\log(c) - T} > 0
		s += ld->c[i] - opt->thres;
		sl = s - ld->f[i];
		if (sl > 0.0) break;
	} // here, i-1+k is the minimum detectable homopolymer length
	ld->max_test = (int)(i * log(i) / opt->thres);
	return ld;
}

void ld_data_destroy(ld_data_t *ld)
{
	void *km = ld->km;
	kfree(km, ld->intv);
	kfree(km, ld->for_pos);
	kdq_destroy(uint32_t, ld->q);
	kfree(km, ld->c); kfree(km, ld->f);
	kfree(km, ld->ht);
	kfree(km, ld);
}

static int32_t ld_dust_forward(ld_data_t *ld, int32_t i0, double max_back, uint16_t *ht)
{
	const ld_opt_t *opt = ld->opt;
	int32_t max_i = -1, i, l;
	double s, sl, max_sf = 0.0;
	memset(ht, 0, sizeof(uint16_t) * (1U<<2*opt->kmer));
	for (i = i0, l = 1, s = sl = 0.0; i < kdq_size(ld->q); ++i, ++l) {
		uint32_t x = kdq_at(ld->q, i);
		s += (x&1? 0 : ld->c[++ht[x>>1]]) - opt->thres;
		sl = s - ld->f[l];
		if (sl >= max_sf) max_sf = sl, max_i = i;
		if (sl > max_back + 1e-6) break;
	}
	return max_i;
}

static int32_t ld_dust_backward(ld_data_t *ld, int64_t pos, const uint16_t *win_ht, double win_sum)
{
	const ld_opt_t *opt = ld->opt;
	double xdrop = opt->thres * (opt->xdrop_len > 0? opt->xdrop_len : opt->ws);
	int32_t i, l, max_i = -1, max_end;
	double s, sl, sw, max_sb = 0.0, last_sl = -1.0;

	memset(ld->ht, 0, sizeof(uint16_t) * (1U<<2*opt->kmer));
	ld->n_for_pos = 0;
	for (i = kdq_size(ld->q) - 1, l = 1, s = sl = sw = 0.0; i >= 0; --i, ++l) { // backward
		uint32_t x = kdq_at(ld->q, i);
		s += (x&1? 0 : ld->c[++ld->ht[x>>1]]) - opt->thres;
		sl = s - ld->f[l]; // this is the score
		sw += (x&1? 0 : ld->c[win_ht[x>>1] + 1 - ld->ht[x>>1]]) - opt->thres;
		if (sw - ld->f[l] < 0.0) break; // in this case, the forward pass won't reach _pos_
		if (sl < last_sl && last_sl > 0.0 && last_sl == max_sb) // store positions where forward may be needed
			ld->for_pos[ld->n_for_pos].pos = i + 1, ld->for_pos[ld->n_for_pos++].max = max_sb;
		if (sl >= max_sb) {
			max_sb = sl, max_i = i;
		} else if (max_i >= 0) {
			if (max_sb - sl > xdrop) break; // X-drop
		}
		//if (win_sum - ld->f[l] - l * opt->thres < max_sb) break; // in this case, we won't get a higher score even if we reach i==0; not very effective
		last_sl = sl;
	}
	if (max_i < 0) return -1;
	if (ld->n_for_pos == 0 || max_i < ld->for_pos[ld->n_for_pos - 1].pos) // this may happen when the max_sb is achieved at the last cycle
		ld->for_pos[ld->n_for_pos].pos = max_i, ld->for_pos[ld->n_for_pos++].max = max_sb;
	for (i = ld->n_for_pos - 1, max_end = -1; i >= 0; --i) { // forward
		const ld_forpos_t *p = &ld->for_pos[i];
		int32_t k;
		if (p->pos < max_end) continue;
		k = ld_dust_forward(ld, p->pos, p->max, ld->ht);
		if (k == kdq_size(ld->q) - 1) return p->pos;
		if (opt->approx) break; // in the approximate mode, do one forward pass only
		max_end = max_end > k? max_end : k;
	}
	return -1;
}

static inline int32_t ld_extend(ld_data_t *ld)
{
	uint32_t x = kdq_last(ld->q);
	int32_t l = kdq_size(ld->q) - 1;
	double diff;
	if (x&1) return -1;
	diff = ld->c[ld->ht[x>>1] + 1] - (ld->f[l + 1] - ld->f[l]);
	if (diff < ld->opt->thres) return -1; // if this doesn't increase score, don't extend
	++ld->ht[x>>1];
	return 0;
}

static int32_t ld_if_backward(const ld_data_t *ld, const uint16_t *ht, int32_t max_step)
{
	int32_t i, j;
	double s = 0.0;
	for (i = kdq_size(ld->q) - 1, j = 0; i >= 0 && j < max_step; --i, ++j) {
		uint32_t x = kdq_at(ld->q, i);
		s += (x&1? 0 : ld->c[ht[x>>1]]) - ld->opt->thres;
		if (s < 0.0) return 0;
	}
	return 1;
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
	int64_t last_q = -1; // last_q: pos in queue when last_i is set
	uint16_t *ht;
	double ht_sum = 0.0;

	ld->n_intv = 0;
	ld->q->front = ld->q->count = 0;
	ht = Kcalloc(ld->km, uint16_t, mask + 1);
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
			if (last_q == 0) {
				if ((*p&1) == 0 && ld->ht[*p>>1] > 0) --ld->ht[*p>>1]; // FIXME: why ld->ht[*p>>1] can be 0?
			} else --last_q; // this needs to be updated as the queue is shifted
		}
		kdq_push(uint32_t, ld->q, x<<1|ambi);
		if (ambi) continue;
		ht_sum += ld->c[++ht[x]];

		j = -1;
		if (ht[x] >= opt->min_start_cnt) { // no need to call the following if x has low count in the window; DON'T test ld_is_back() here!
			double swin = ht_sum - ld->f[kdq_size(ld->q)] - kdq_size(ld->q) * opt->thres; // this is the full window score
			if (i == en && (last_q == 0 || i - st >= kdq_size(ld->q)) && swin > 0.0) // test and potentially extend the base at i, ONLY when the full window is LCR
				j = ld_extend(ld);
			if (j < 0 && ld_if_backward(ld, ht, ld->max_test))
				j = ld_dust_backward(ld, i, ht, ht_sum);
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
			last_q = j;
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
