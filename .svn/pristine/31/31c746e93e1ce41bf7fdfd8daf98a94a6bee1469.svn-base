#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Given a distance matrix, computes connectivity length.
 *
 * Marchiori and Latora (2000). Harmony in the small-world. Physica A 285:
 * 539-546.
 */
FP_T BCT_NAMESPACE::connectivity_length(const MATRIX_T* D) {
	if (safe_mode) check_status(D, SQUARE, "connectivity_length");
	int N = D->size1;
	FP_T sum = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				continue;
			}
			FP_T value = MATRIX_ID(get)(D, i, j);
			if (gsl_finite(value) == 1 && fp_nonzero(value)) {
				sum += 1.0 / value;
			}
		}
	}
	if (fp_zero(sum)) {
		return GSL_POSINF;
	} else {
		return (FP_T)(N * (N - 1)) / sum;
	}
}
