#include "bct.h"

/*
 * Computes the clustering coefficient for a binary undirected graph.
 */
VECTOR_T* BCT_NAMESPACE::clustering_coef_bu(const MATRIX_T* G) {
	if (safe_mode) check_status(G, SQUARE | BINARY | UNDIRECTED, "clustering_coef_bu");
	
	// n=length(G);
	int n = length(G);
	
	// C=zeros(n,1);
	VECTOR_T* C = zeros_vector(n);
	
	// for u=1:n
#ifdef _OPENMP
#pragma omp parallel for shared(C)
#endif
	for (int u = 0; u < n; u++) {
		
		// V=find(G(u,:));
		// k=length(V);
		// if k>=2;
		VECTOR_ID(const_view) G_row_u = MATRIX_ID(const_row)(G, u);
		VECTOR_T* V = find(&G_row_u.vector);
		if (V != NULL) {
			int k = length(V);
			if (k >= 2) {
				
				// S=G(V,V);
				MATRIX_T* S = ordinal_index(G, V, V);
				
				// C(u)=sum(S(:))/(k^2-k);
				VECTOR_T* sum_S = sum(S);
				MATRIX_ID(free)(S);
				VECTOR_ID(set)(C, u, sum(sum_S) / (FP_T)(k * (k - 1)));
				VECTOR_ID(free)(sum_S);
			}
			
			VECTOR_ID(free)(V);
		}
	}
	
	return C;
}
