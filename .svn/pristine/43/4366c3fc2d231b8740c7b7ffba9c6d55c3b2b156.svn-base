#include "bct.h"

/*
 * Computes strength, in-strength, and out-strength for a directed graph.
 */
VECTOR_T* BCT_NAMESPACE::strengths_dir(const MATRIX_T* CIJ, VECTOR_T** is, VECTOR_T** os) {
	if (safe_mode) check_status(CIJ, SQUARE | DIRECTED, "strengths_dir");
	
	// is = sum(CIJ,1);
	VECTOR_T* _is = sum(CIJ, 1);
	
	// os = sum(CIJ,2);
	VECTOR_T* _os = sum(CIJ, 2);
	
	// str = is+os;
	VECTOR_T* str = copy(_is);
	VECTOR_ID(add)(str, _os);
	
	if (is != NULL) *is = _is; else VECTOR_ID(free)(_is);
	if (os != NULL) *os = _os; else VECTOR_ID(free)(_os);
	return str;
}
