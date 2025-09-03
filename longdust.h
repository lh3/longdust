#ifndef LONGDUST_H
#define LONGDUST_H

#define LD_VERSION "1.0-r51-dirty"

#include <stdint.h>

typedef struct {
	int32_t kmer, ws, xdrop_len, exact;
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

/**
 * Initialize default parameters
 */
void ld_opt_init(ld_opt_t *opt);

/**
 * Initialize the internal buffer for LCR identification
 *
 * @param km       NULL or memory space returned by km_init()
 * @param opt      setting
 *
 * @return initialized buffer
 */
ld_data_t *ld_data_init(void *km, const ld_opt_t *opt);

/**
 * Deallocate the buffer
 */
void ld_data_destroy(ld_data_t *ld);

/**
 * Find LCRs on the forward strand
 */
void ld_dust1(ld_data_t *ld, int64_t len, const uint8_t *seq);

/**
 * The union of LCRs found on both strands
 */
void ld_dust2(ld_data_t *ld, int64_t len, const uint8_t *seq);

/**
 * Get the list of LCRs stored in the buffer
 *
 * @param ld        result buffer
 * @param n         number of intervals (out)
 *
 * @return array of intervals of length _n_ (read-only)
 */
const ld_intv_t *ld_get_intv(const ld_data_t *ld, int64_t *n);

#ifdef __cplusplus
}
#endif

#endif
