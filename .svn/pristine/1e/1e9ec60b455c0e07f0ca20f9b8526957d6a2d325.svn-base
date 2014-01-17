#include <gsl/gsl_math.h>

#include "bct.h"

MATRIX_T* distance_inv(const MATRIX_T*, const MATRIX_T*);

/*
 * Computes global efficiency.  Takes an optional distance matrix that is
 * computed if not given.
 */
FP_T BCT_NAMESPACE::efficiency_global(const MATRIX_T* G, const MATRIX_T* D) {
	if (safe_mode) check_status(G, SQUARE, "efficiency_global");
	
	// N=length(G);
	int N = length(G);
	
	// e=distance_inv(G);
	MATRIX_T* e = distance_inv(G, D);
	
	// E=sum(e(:))./(N^2-N);
	VECTOR_T* e_v = to_vector(e);
	MATRIX_ID(free)(e);
	FP_T sum_e = sum(e_v);
	VECTOR_ID(free)(e_v);
	return sum_e / (FP_T)(N * (N - 1));
}

/*
 * Computes local efficiency.
 */
VECTOR_T* BCT_NAMESPACE::efficiency_local(const MATRIX_T* G) {
	if (safe_mode) check_status(G, SQUARE, "efficiency_local");
	
	// N=length(G);
	int N = length(G);
	
	// E=zeros(N,1);
	VECTOR_T* E = zeros_vector(N);
	
	// for u=1:N
#ifdef _OPENMP
#pragma omp parallel for shared(E)
#endif
	for (int u = 0; u < N; u++) {
		
		// V=find(G(u,:));
		VECTOR_ID(const_view) G_row_u = MATRIX_ID(const_row)(G, u);
		VECTOR_T* V = find(&G_row_u.vector);
		if (V != NULL) {
			
			// k=length(V);
			int k = length(V);
			
			// if k>=2;
			if (k >= 2) {
				
				// e=distance_inv(G(V,V));
				MATRIX_T* G_idx = ordinal_index(G, V, V);
				MATRIX_T* e = distance_inv(G_idx, NULL);
				MATRIX_ID(free)(G_idx);
				
				// E(u)=sum(e(:))./(k^2-k);
				VECTOR_T* e_v = to_vector(e);
				MATRIX_ID(free)(e);
				FP_T sum_e = sum(e_v);
				VECTOR_ID(free)(e_v);
				VECTOR_ID(set)(E, u, sum_e / (FP_T)(k * (k - 1)));
			}
			
			VECTOR_ID(free)(V);
		}
	}
	
	return E;
}

MATRIX_T* distance_inv(const MATRIX_T* G, const MATRIX_T* D) {
	using namespace BCT_NAMESPACE;
	
	MATRIX_T* D_inv;
	if (D == NULL) {
		MATRIX_T* G_inv = invert_elements(G);
		MATRIX_T* temp = distance_wei(G_inv);
		MATRIX_ID(free)(G_inv);
		D_inv = invert_elements(temp);
		MATRIX_ID(free)(temp);
	} else {
		D_inv = invert_elements(D);
	}
	return D_inv;
}
