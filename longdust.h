#ifndef LONGDUST_H
#define LONGDUST_H

#include <stdint.h>

typedef struct {
	int32_t kmer, ws;
	double thres, xdrop;
} ld_opt_t;

typedef struct {
	int64_t st, en;
} ld_intv_t;

struct ld_data_s;
typedef struct ld_data_s ld_data_t;

#endif
