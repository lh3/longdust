#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "kalloc.h"
#include "khashl.h"
#include "kdq.h"

double ld_f_large(double lambda)
{
	double x = 0.5 * log(2.0 * M_PI * M_E * lambda) - 1.0 / 12.0 / lambda * (1.0 + 0.5 / lambda + 19.0 / 30.0 / lambda / lambda);
	x += lambda * (log(lambda) - 1.0);
	return x;
}

double *ld_cal_f(void *km, int32_t k, int32_t max_l, double eps)
{
	static const int32_t max_n = 1000;
	double *f;
	int32_t l;
	assert(k < 16);
	f = Kcalloc(km, double, max_l + 1);
	for (l = 1; l <= max_l; ++l) {
		double lambda = (double)l / (1U<<2*k);
		double x = 0.0, sn = 0.0, y = lambda;
		int32_t n;
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
	}
	return f;
}

int main(int argc, char *argv[])
{
	double *f;
	int i;
	f = ld_cal_f(0, 3, 1000, 1e-6);
//	for (i = 0; i < 1000; ++i)
//		printf("%d\t%.9f\n", i, f[i]);
	return 0;
}
