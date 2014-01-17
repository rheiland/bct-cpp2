#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * WARNING: BCT_NAMESPACE::charpath_lambda takes a distance matrix, but
 * BCT_NAMESPACE::capped_charpath_lambda takes a connection matrix.  Both should be
 * lengths, not weights (called distances in CalcMetric).
 */

/*
 * Given a distance matrix, computes characteristic path length.
 */
FP_T BCT_NAMESPACE::charpath_lambda(const MATRIX_T* D) {
	if (safe_mode) check_status(D, SQUARE, "charpath_lambda");
	
	// lambda = sum(sum(D(D~=Inf)))/length(nonzeros(D~=Inf));
	MATRIX_T* D_neq_inf = compare_elements(D, fp_not_equal, GSL_POSINF);
	VECTOR_T* D_idx = logical_index(D, D_neq_inf);
	FP_T sum_D_idx = sum(D_idx);
	VECTOR_ID(free)(D_idx);
	FP_T ret = sum_D_idx / (FP_T)nnz(D_neq_inf);
	MATRIX_ID(free)(D_neq_inf);
	return ret;
}

/*
 * Given a connection matrix, computes capped characteristic path length.
 */
FP_T BCT_NAMESPACE::capped_charpath_lambda(const MATRIX_T* L) {
	if (safe_mode) check_status(L, SQUARE, "capped_charpath_lambda");
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
	FP_T dmax = (FP_T)N * lmean;
	FP_T dmean = 0.0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				continue;
			}
			FP_T d = MATRIX_ID(get)(D, i, j);
			dmean += (d < dmax) ? d : dmax;
		}
	}
	dmean /= N * (N - 1);
	MATRIX_ID(free)(D);
	return dmean;
}

/*
 * Given a distance matrix, computes eccentricity, radius, and diameter.
 */
VECTOR_T* BCT_NAMESPACE::charpath_ecc(const MATRIX_T* D, FP_T* radius, FP_T* diameter) {
	if (safe_mode) check_status(D, SQUARE, "charpath_ecc");
	
	// ecc = max(D.*(D~=Inf),[],2);
	MATRIX_T* D_finite = copy(D);
	MATRIX_T* D_eq_inf = compare_elements(D, fp_equal, GSL_POSINF);
	logical_index_assign(D_finite, D_eq_inf, 0.0);
	MATRIX_ID(free)(D_eq_inf);
	VECTOR_T* ecc = max(D_finite, 2);
	MATRIX_ID(free)(D_finite);
	
	// radius = min(ecc);
	if (radius != NULL) {
		*radius = min(ecc);
	}
	
	// diameter = max(ecc);
	if (diameter != NULL) {
		*diameter = max(ecc);
	}
	
	return ecc;
}
