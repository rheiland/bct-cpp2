#include "bct.h"

/*
 * Computes density for an undirected graph.  Connection weights are ignored.
 */
FP_T BCT_NAMESPACE::density_und(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE | UNDIRECTED, "density_und");
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// K = nnz(triu(CIJ));
	MATRIX_T* triu_CIJ = triu(CIJ);
	int K = nnz(triu_CIJ);
	MATRIX_ID(free)(triu_CIJ);
	
	// kden = K/((N^2-N)/2);
	return (FP_T)K / ((FP_T)(N * (N - 1)) / 2.0);
}
