#include "bct.h"

/*
 * Computes degree for an undirected graph.  Connection weights are ignored.
 */
VECTOR_T* BCT_NAMESPACE::degrees_und(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE | UNDIRECTED, "degrees_und");
	
	// CIJ = double(CIJ~=0);
	// deg = sum(CIJ);
	VECTOR_T* deg = VECTOR_ID(alloc)(CIJ->size2);

#ifdef _OPENMP
#pragma omp parallel for shared(deg)
#endif
	for (int i = 0; i < (int)CIJ->size2; i++) {
		VECTOR_ID(const_view) CIJ_col_i = MATRIX_ID(const_column)(CIJ, i);
		VECTOR_ID(set)(deg, i, nnz(&CIJ_col_i.vector));
	}
	return deg;
}
