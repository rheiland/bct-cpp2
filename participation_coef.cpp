#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Computes nodal participation coefficient and community structure.  For a
 * directed graph, computes "out-neighbor" participation coefficient.
 */
VECTOR_T* BCT_NAMESPACE::participation_coef(const MATRIX_T* W, const VECTOR_T* Ci) {
	if (safe_mode) check_status(W, SQUARE | BINARY, "participation_coef");
	
	// n=length(W);
	int n = length(W);
	
	// Ko=sum(W,2);
	VECTOR_T* Ko = sum(W, 2);
	
	// Gc=(W~=0)*diag(Ci);
	MATRIX_T* W_neq_0 = compare_elements(W, fp_not_equal, 0.0);
	MATRIX_T* diag_Ci = diag(Ci);
	MATRIX_T* Gc = mul(W_neq_0, diag_Ci);
	MATRIX_ID(free)(W_neq_0);
	MATRIX_ID(free)(diag_Ci);
	
	// Kc2=zeros(n,1);
	VECTOR_T* Kc2 = zeros_vector(n);
	
	// for i=1:max(Ci);
	for (int i = 1; i <= (int)max(Ci); i++) {
		
		// Kc2=Kc2+(sum(W.*(Gc==i),2).^2);
		MATRIX_T* Gc_eq_i = compare_elements(Gc, fp_equal, (FP_T)i);
		MATRIX_T* W_mul_Gc_eq_i = copy(W);
		MATRIX_ID(mul_elements)(W_mul_Gc_eq_i, Gc_eq_i);
		MATRIX_ID(free)(Gc_eq_i);
		VECTOR_T* sum_W_mul_Gc_eq_i = sum(W_mul_Gc_eq_i, 2);
		MATRIX_ID(free)(W_mul_Gc_eq_i);
		VECTOR_T* sum_W_mul_Gc_eq_i_pow_2 = pow_elements(sum_W_mul_Gc_eq_i, 2.0);
		VECTOR_ID(free)(sum_W_mul_Gc_eq_i);
		VECTOR_ID(add)(Kc2, sum_W_mul_Gc_eq_i_pow_2);
		VECTOR_ID(free)(sum_W_mul_Gc_eq_i_pow_2);
	}
	
	MATRIX_ID(free)(Gc);
	
	// P=ones(n,1)-Kc2./(Ko.^2);
	VECTOR_T* Ko_pow_2 = pow_elements(Ko, 2.0);
	VECTOR_ID(div)(Kc2, Ko_pow_2);
	VECTOR_ID(free)(Ko_pow_2);
	VECTOR_T* P = ones_vector(n);
	VECTOR_ID(sub)(P, Kc2);
	VECTOR_ID(free)(Kc2);
	
	// P(~Ko)=0;
	VECTOR_T* not_Ko = logical_not(Ko);
	VECTOR_ID(free)(Ko);
	logical_index_assign(P, not_Ko, 0.0);
	VECTOR_ID(free)(not_Ko);
	
	return P;
}
