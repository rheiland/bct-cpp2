#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Computes the distance matrix for a binary graph.
 */
MATRIX_T* BCT_NAMESPACE::distance_bin(const MATRIX_T* G) {
	if (safe_mode) check_status(G, SQUARE | BINARY, "distance_bin");
	
	// D=eye(length(G));
	MATRIX_T* D = eye(length(G));
	
	// n=1;
	int n = 1;
	
	// nPATH=G;
	MATRIX_T* nPATH = copy(G);
	
	// L=(nPATH~=0);
	MATRIX_T* L = compare_elements(nPATH, fp_not_equal, 0.0);
	
	// while find(L,1);
	VECTOR_T* find_L = find(L, 1);
	while (find_L != NULL) {
		VECTOR_ID(free)(find_L);
		
		// D=D+n.*L;
		MATRIX_ID(scale)(L, (FP_T)n);
		MATRIX_ID(add)(D, L);
		
		// n=n+1;
		n++;
		
		// nPATH=nPATH*G;
		MATRIX_T* temp = mul(nPATH, G);
		MATRIX_ID(free)(nPATH);
		nPATH = temp;
		
		// L=(nPATH~=0).*(D==0);
		MATRIX_ID(free)(L);
		L = compare_elements(nPATH, fp_not_equal, 0.0);
		MATRIX_T* D_eq_0 = compare_elements(D, fp_equal, 0.0);
		MATRIX_ID(mul_elements)(L, D_eq_0);
		MATRIX_ID(free)(D_eq_0);
		
		find_L = find(L, 1);
	}
	
	MATRIX_ID(free)(nPATH);
	MATRIX_ID(free)(L);
	
	// D(~D)=inf;
	MATRIX_T* not_D = logical_not(D);
	logical_index_assign(D, not_D, GSL_POSINF);
	MATRIX_ID(free)(not_D);
	
	// D=D-eye(length(G));
	MATRIX_T* eye_length_G = eye(length(G));
	MATRIX_ID(sub)(D, eye_length_G);
	MATRIX_ID(free)(eye_length_G);
	
	return D;
}
