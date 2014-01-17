#include "bct.h"

/*
 * Computes density for a directed graph.  Connection weights are ignored.
 */
FP_T BCT_NAMESPACE::density_dir(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE | DIRECTED, "density_dir");
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// K = nnz(CIJ);
	int K = nnz(CIJ);
	
	// kden = K/(N^2-N);
	return (FP_T)K / (FP_T)(N * (N - 1));
}
