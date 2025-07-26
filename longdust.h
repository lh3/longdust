#ifndef LONGDUST_H
#define LONGDUST_H

#include <stdint.h>

typedef struct {
	int32_t kmer, ws;
	double thres, xdrop;
} ld_opt_t;

#endif
