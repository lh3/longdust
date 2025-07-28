#ifndef LONGDUST_H
#define LONGDUST_H

#define LD_VERSION "r24"

#include <stdint.h>

typedef struct {
	int32_t kmer, ws, xdrop_len1, xdrop_len2;
	double thres;
} ld_opt_t;

typedef struct {
	int64_t st, en;
} ld_intv_t;

struct ld_data_s;
typedef struct ld_data_s ld_data_t;

#ifdef __cplusplus
extern "C" {
#endif

void ld_opt_init(ld_opt_t *opt);
ld_data_t *ld_data_init(void *km, const ld_opt_t *opt);
void ld_data_destroy(ld_data_t *ld);

void ld_dust1(ld_data_t *ld, int64_t len, const uint8_t *seq);
void ld_dust2(ld_data_t *ld, int64_t len, const uint8_t *seq);
const ld_intv_t *ld_get_intv(const ld_data_t *ld, int64_t *n);

#ifdef __cplusplus
}
#endif

#endif
