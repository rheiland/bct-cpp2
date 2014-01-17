#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Computes reachability and distance matrices using breadth-first search.
 */
MATRIX_T* BCT_NAMESPACE::breadthdist(const MATRIX_T* CIJ, MATRIX_T** D) {
	if (safe_mode) check_status(CIJ, SQUARE, "breadthdist");
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// D = zeros(N);
	MATRIX_T* _D = zeros(N);
	
	// for i=1:N
	for (int i = 0; i < N; i++) {
		
		// D(i,:) = breadth(CIJ,i);
		VECTOR_T* distance = breadth(CIJ, i);
		MATRIX_ID(set_row)(_D, i, distance);
		VECTOR_ID(free)(distance);
	}
	
	// D(D==0) = Inf;
	MATRIX_T* D_eq_0 = compare_elements(_D, fp_equal, 0.0);
	logical_index_assign(_D, D_eq_0, GSL_POSINF);
	MATRIX_ID(free)(D_eq_0);
	
	// R = double(D~=Inf);
	MATRIX_T* R = compare_elements(_D, fp_not_equal, GSL_POSINF);
	
	if (D != NULL) *D = _D; else MATRIX_ID(free)(_D);
	return R;
}
