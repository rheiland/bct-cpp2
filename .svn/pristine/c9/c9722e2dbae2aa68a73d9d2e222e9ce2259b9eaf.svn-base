#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Computes clustering coefficient for a binary directed graph.
 */
VECTOR_T* BCT_NAMESPACE::clustering_coef_bd(const MATRIX_T* A) {
	if (safe_mode) check_status(A, SQUARE | BINARY | DIRECTED, "clustering_coef_bd");
	
	// S=A+A.';
	MATRIX_T* S = MATRIX_ID(alloc)(A->size2, A->size1);
	MATRIX_ID(transpose_memcpy)(S, A);
	MATRIX_ID(add)(S, A);
	
	// K=sum(S,2);
	VECTOR_T* K = sum(S, 2);
	
	// cyc3=diag(S^3)/2;
	MATRIX_T* S_pow_3_div_2 = pow(S, 3);
	MATRIX_ID(free)(S); 
	MATRIX_ID(scale)(S_pow_3_div_2, 0.5);
	VECTOR_ID(view) cyc3 = MATRIX_ID(diagonal)(S_pow_3_div_2);
	
	// K(cyc3==0)=inf;
	VECTOR_T* cyc3_eq_0 = compare_elements(&cyc3.vector, fp_equal, 0.0);
	logical_index_assign(K, cyc3_eq_0, GSL_POSINF);
	VECTOR_ID(free)(cyc3_eq_0);
	
	// CYC3=K.*(K-1)-2*diag(A^2);
	VECTOR_T* K_sub_1 = copy(K);
	VECTOR_ID(add_constant)(K_sub_1, -1.0);
	MATRIX_T* A_pow_2_mul_2 = pow(A, 2);
	MATRIX_ID(scale)(A_pow_2_mul_2, 2.0);
	VECTOR_ID(view) diag_A_pow_2_mul_2 = MATRIX_ID(diagonal)(A_pow_2_mul_2);
	VECTOR_T* CYC3 = K;
	VECTOR_ID(mul)(CYC3, K_sub_1);
	VECTOR_ID(free)(K_sub_1);
	VECTOR_ID(sub)(CYC3, &diag_A_pow_2_mul_2.vector);
	MATRIX_ID(free)(A_pow_2_mul_2);
	
	// C=cyc3./CYC3
	VECTOR_T* C = copy(&cyc3.vector);
	MATRIX_ID(free)(S_pow_3_div_2);
	VECTOR_ID(div)(C, CYC3);
	VECTOR_ID(free)(CYC3);
	
	return C;
}
