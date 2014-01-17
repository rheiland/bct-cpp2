#include <cmath>

#include "bct.h"

/*
 * Finds the dominant eigenvector using power iteration.  Adapted from
 * NetworkX v1.4 eigenvector_centrality.
 *
 * http://www.mathworks.de/matlabcentral/fx_files/7978/1/mPowerEig.c
 * http://en.wikipedia.org/wiki/Power_iteration
 */
VECTOR_T* BCT_NAMESPACE::eigenvector_centrality(const MATRIX_T* G) {
	if (safe_mode) check_status(G, SQUARE, "eigenvector_centrality");
	
	FP_T tol = 1e-6;  // Ensures that the average value is stable at 1e-6
	int maxiter = 1000;
	int N = G->size1;
	VECTOR_T* x = VECTOR_ID(alloc)(N);
	VECTOR_T* xlast = VECTOR_ID(alloc)(N);
	
	// Note that this starting vector is already normalized (elements sum to 1)
	FP_T startval = 1.0 / (FP_T)N;
	VECTOR_ID(set_all)(x, startval);
	
	FP_T evec_norm;
	FP_T err;
	FP_T newval;
	for (int iter = 0; iter < maxiter; iter++) {
		err = 0.0;
		VECTOR_ID(memcpy)(xlast, x);
		VECTOR_ID(set_zero)(x);
		
#ifdef _OPENMP
#pragma omp parallel for private(newval) shared(x, xlast)
#endif
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				newval = VECTOR_ID(get)(x, i) + VECTOR_ID(get)(xlast, j) * MATRIX_ID(get)(G, i, j);
				VECTOR_ID(set)(x, i, newval);
			}
		}
		
		// Normalize vector
		evec_norm = norm(x, 2);
		VECTOR_ID(scale)(x, 1.0 / evec_norm);
		
		// Check convergence
		for (int i = 0; i < (int)x->size; i++) {
			err += std::abs(VECTOR_ID(get)(x, i) - VECTOR_ID(get)(xlast, i));
		}
		if (err < (FP_T)N * tol) {
			break;  // End power iteration
		}
		
	}
	VECTOR_ID(free)(xlast);
	return x;
}
