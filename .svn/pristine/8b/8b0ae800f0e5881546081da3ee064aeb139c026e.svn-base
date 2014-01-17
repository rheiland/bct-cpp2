#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Computes z-score for a binary graph and its corresponding community
 * structure.  For a directed graph, computes out-neighbor z-score.
 */
VECTOR_T* BCT_NAMESPACE::module_degree_zscore(const MATRIX_T* A, const VECTOR_T* Ci) {
	if (safe_mode) check_status(A, SQUARE | BINARY, "module_degree_zscore");
	
	// n=length(A);
	int n = length(A);
	
	// Z=zeros(n,1);
	VECTOR_T* Z = zeros_vector(n);
	
	// for i=1:max(Ci)
	for (int i = 1; i <= (int)max(Ci); i++) {
		
		// Koi=sum(A(Ci==i,Ci==i),2);
		VECTOR_T* Ci_eq_i = compare_elements(Ci, fp_equal, (FP_T)i);
		MATRIX_T* A_idx = logical_index(A, Ci_eq_i, Ci_eq_i);
		VECTOR_T* Koi = sum(A_idx, 2);
		MATRIX_ID(free)(A_idx);
		
		// Z(Ci==i)=(Koi-mean(Koi))./std(Koi);
		FP_T std_Koi = MATLAB_NAMESPACE::std(Koi);
		VECTOR_ID(add_constant)(Koi, -mean(Koi));
		VECTOR_ID(scale)(Koi, 1.0 / std_Koi);
		logical_index_assign(Z, Ci_eq_i, Koi);
		VECTOR_ID(free)(Ci_eq_i);
		VECTOR_ID(free)(Koi);
	}
	
	// Z(isnan(Z))=0;
	for (int i = 0; i < (int)Z->size; i++) {
		if (gsl_isnan(VECTOR_ID(get)(Z, i)) == 1) {
			VECTOR_ID(set)(Z, i, 0.0);
		}
	}
	
	return Z;
}
