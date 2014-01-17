#include "bct.h"

/*
 * Computes strength for an undirected graph.
 */
VECTOR_T* BCT_NAMESPACE::strengths_und(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE | UNDIRECTED, "strengths_und");
	
	// str = sum(CIJ);
	return sum(CIJ);
}
