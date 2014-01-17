#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Computes the clustering coefficient for a weighted undirected graph.
 */
VECTOR_T* BCT_NAMESPACE::clustering_coef_wu(const MATRIX_T* W) {
	if (safe_mode) check_status(W, SQUARE | WEIGHTED | UNDIRECTED, "clustering_coef_wu");
	
	// K=sum(W~=0,2);
	MATRIX_T* W_neq_0 = compare_elements(W, fp_not_equal, 0.0);
	VECTOR_T* K = sum(W_neq_0, 2);
	MATRIX_ID(free)(W_neq_0);
	
	// cyc3=diag((W.^(1/3))^3);
	MATRIX_T* W_pow_1_3 = pow_elements(W, 1.0 / 3.0);
	MATRIX_T* W_pow_1_3_pow_3 = pow(W_pow_1_3, 3);
	MATRIX_ID(free)(W_pow_1_3);
	VECTOR_ID(view) cyc3 = MATRIX_ID(diagonal)(W_pow_1_3_pow_3);
	
	// K(cyc3==0)=inf;
	VECTOR_T* cyc3_eq_0 = compare_elements(&cyc3.vector, fp_equal, 0.0);
	logical_index_assign(K, cyc3_eq_0, GSL_POSINF);
	VECTOR_ID(free)(cyc3_eq_0);
	
	// C=cyc3./(K.*(K-1));
	VECTOR_T* K_sub_1 = copy(K);
	VECTOR_ID(add_constant)(K_sub_1, -1.0);
	VECTOR_ID(mul)(K, K_sub_1);
	VECTOR_ID(free)(K_sub_1);
	VECTOR_T* C = copy(&cyc3.vector);
	MATRIX_ID(free)(W_pow_1_3_pow_3);
	VECTOR_ID(div)(C, K);
	VECTOR_ID(free)(K);
	
	return C;
}
