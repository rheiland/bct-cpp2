#include "bct.h"

/*
 * Generates a random undirected binary graph with N nodes and K edges.  No
 * edges are placed on the main diagonal.
 */
MATRIX_T* BCT_NAMESPACE::makerandCIJ_bu(int N, int K) {
	
	// ind = triu(~eye(N));
	MATRIX_T* eye_N = eye(N);
	MATRIX_T* not_eye_N = logical_not(eye_N);
	MATRIX_ID(free)(eye_N);
	MATRIX_T* ind = triu(not_eye_N);
	MATRIX_ID(free)(not_eye_N);
	
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
	
	// CIJ = CIJ+CIJ';
	MATRIX_T* CIJ_transpose = MATRIX_ID(alloc)(CIJ->size2, CIJ->size1);
	MATRIX_ID(transpose_memcpy)(CIJ_transpose, CIJ);
	MATRIX_ID(add)(CIJ, CIJ_transpose);
	MATRIX_ID(free)(CIJ_transpose);
	
	return CIJ;
}
