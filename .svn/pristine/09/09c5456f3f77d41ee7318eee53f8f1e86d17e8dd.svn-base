#include <cmath>

#include "bct.h"

/*
 * WARNING: BCT_NAMESPACE::normalized_path_length takes a distance matrix, but
 * BCT_NAMESPACE::normalized_path_length_m takes a connection matrix.  Both should be
 * lengths, not weights (called distances in CalcMetric).
 */

/*
 * Given a distance matrix, computes the normalized shortest path length.
 */
FP_T BCT_NAMESPACE::normalized_path_length(const MATRIX_T* D, FP_T wmax) {
	if (safe_mode) check_status(D, SQUARE, "normalized_path_length");
	int N = D->size1;
	FP_T dmin = 1.0 / wmax;
	FP_T dmax = (FP_T)N / wmax;
	FP_T sum = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				continue;
			}
			FP_T d = MATRIX_ID(get)(D, i, j);
			sum += (d < dmax) ? d : dmax;
		}
	}
	return std::abs(((sum / (FP_T)(N * (N - 1))) - dmin) / (dmax - dmin));
}

/*
 * Computes the normalized shortest path length using dmax = N * lmean, where
 * lmean is the average distance between all directly connected nodes.
 */
FP_T BCT_NAMESPACE::normalized_path_length_m(const MATRIX_T* L, FP_T wmax) {
	if (safe_mode) check_status(L, SQUARE, "normalized_path_length_m");
	int N = L->size1;
	int nonzeros = 0;
	FP_T lmean = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				continue;
			}
			FP_T l = MATRIX_ID(get)(L, i, j);
			if (fp_nonzero(l)) {
				nonzeros++;
				lmean += l;
			}
		}
	}
	lmean /= nonzeros;
	MATRIX_T* D = distance_wei(L);
	FP_T dmin = 1.0 / wmax;
	FP_T dmax = (FP_T)N * lmean;
	FP_T sum = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				continue;
			}
			FP_T d = MATRIX_ID(get)(D, i, j);
			sum += (d < dmax) ? d : dmax;
		}
	}
	MATRIX_ID(free)(D);
	return std::abs(((sum / (FP_T)(N * (N - 1))) - dmin) / (dmax - dmin));
}
