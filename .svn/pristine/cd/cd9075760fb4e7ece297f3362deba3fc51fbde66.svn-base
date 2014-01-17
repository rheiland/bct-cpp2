#include "bct.h"

/*
 * Applies an absolute weight threshold to a graph.  All weights below this
 * threshold, as well as those on the main diagonal, are set to zero.
 */
MATRIX_T* BCT_NAMESPACE::threshold_absolute(const MATRIX_T* W, FP_T thr) {
	if (safe_mode) check_status(W, SQUARE, "threshold_absolute");
	
	MATRIX_T* W_thr = copy(W);
	
	// W(1:size(W,1)+1:end)=0;
	for (int i = 0; i < (int)W_thr->size1; i++) {
		MATRIX_ID(set)(W_thr, i, i, 0.0);
	}
	
	// W(W<thr)=0;
	MATRIX_T* W_thr_lt_thr = compare_elements(W_thr, fp_less, thr);
	logical_index_assign(W_thr, W_thr_lt_thr, 0.0);
	MATRIX_ID(free)(W_thr_lt_thr);
	
	return W_thr;
}
