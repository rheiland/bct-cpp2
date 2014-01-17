#include "bct.h"

/*
 * Generates a random directed binary graph with N nodes and K edges.  No edges
 * are placed on the main diagonal.
 */
MATRIX_T* BCT_NAMESPACE::makerandCIJ_bd(int N, int K) {
	
	// ind = ~eye(N);
	MATRIX_T* eye_N = eye(N);
	MATRIX_T* ind = logical_not(eye_N);
	MATRIX_ID(free)(eye_N);
	
	// i = find(ind);
	VECTOR_T* i = find(ind);
	MATRIX_ID(free)(ind);
	
	// rp = randperm(length(i));
	gsl_permutation* rp = randperm(length(i));
	
	// irp = i(rp);
	VECTOR_T* irp = permute(rp, i);
	gsl_permutation_free(rp);
	VECTOR_ID(free)(i);
	
	// CIJ = zeros(N);
	MATRIX_T* CIJ = zeros(N);
	
	// CIJ(irp(1:K)) = 1;
	VECTOR_ID(view) irp_subv = VECTOR_ID(subvector)(irp, 0, K);
	ordinal_index_assign(CIJ, &irp_subv.vector, 1.0);
	VECTOR_ID(free)(irp);
	
	return CIJ;
}
