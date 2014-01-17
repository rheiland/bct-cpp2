#include "bct.h"

/*
 * Computes degree, in-degree, and out-degree for a directed graph.  Connection
 * weights are ignored.
 */
VECTOR_T* BCT_NAMESPACE::degrees_dir(const MATRIX_T* CIJ, VECTOR_T** id, VECTOR_T** od) {
	if (safe_mode) check_status(CIJ, SQUARE | DIRECTED, "degrees_dir");
	
	// CIJ = double(CIJ~=0);
	// id = sum(CIJ,1);
	VECTOR_T* _id = VECTOR_ID(alloc)(CIJ->size2);
	for (int i = 0; i < (int)CIJ->size2; i++) {
		VECTOR_ID(const_view) CIJ_col_i = MATRIX_ID(const_column)(CIJ, i);
		VECTOR_ID(set)(_id, i, nnz(&CIJ_col_i.vector));
	}
	
	// od = sum(CIJ,2);
	VECTOR_T* _od = VECTOR_ID(alloc)(CIJ->size1);
	for (int i = 0; i < (int)CIJ->size1; i++) {
		VECTOR_ID(const_view) CIJ_col_i = MATRIX_ID(const_row)(CIJ, i);
		VECTOR_ID(set)(_od, i, nnz(&CIJ_col_i.vector));
	}
	
	// deg = id+od;
	VECTOR_T* deg = copy(_id);
	VECTOR_ID(add)(deg, _od);
	
	if (id != NULL) *id = _id; else VECTOR_ID(free)(_id);
	if (od != NULL) *od = _od; else VECTOR_ID(free)(_od);
	return deg;
}
