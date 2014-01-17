#include <gsl/gsl_randist.h>

#include "bct.h"

/*
 * Generates a random directed weighted graph with N nodes and K edges.  Weights
 * are chosen uniformly between wmin and wmax.  No edges are placed on the main
 * diagonal.
 */
MATRIX_T* BCT_NAMESPACE::makerandCIJ_wd(int N, int K, FP_T wmin, FP_T wmax) {
	gsl_rng* rng = get_rng();
	VECTOR_T* w = VECTOR_ID(alloc)(K);
	for (int i = 0; i < K; i++) {
		VECTOR_ID(set)(w, i, gsl_rng_uniform(rng) * (wmax - wmin) + wmin);
	}
	
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
	
	// CIJ(irp(1:K)) = w;
	VECTOR_ID(view) irp_subv = VECTOR_ID(subvector)(irp, 0, K);
	ordinal_index_assign(CIJ, &irp_subv.vector, w);
	VECTOR_ID(free)(w);
	VECTOR_ID(free)(irp);
	
	return CIJ;
}

/*
 * Generates a random directed weighted graph with the same number of nodes,
 * number of edges, and weight distribution as the given graph.  No edges are
 * placed on the main diagonal.  The given matrix should therefore not contain
 * nonzero entries on the main diagonal.
 */
MATRIX_T* BCT_NAMESPACE::makerandCIJ_wd_wp(const MATRIX_T* m) {
	if (safe_mode) check_status(m, SQUARE | NO_LOOPS, "makerandCIJ_wd_wp");
	int N = m->size1;
	int K = N * (N - 1);
	FP_T* w = new FP_T[K];
	for (int i = 0, k = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			if (i != j) {
				w[k++] = MATRIX_ID(get)(m, i, j);
			}
		}
	}
	gsl_rng* rng = get_rng();
	gsl_ran_shuffle(rng, w, K, sizeof(FP_T));
	MATRIX_T* rand_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0, k = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			if (i == j) {
				MATRIX_ID(set)(rand_m, i, j, 0.0);
			} else {
				MATRIX_ID(set)(rand_m, i, j, w[k++]);
			}
		}
	}
	delete[] w;
	return rand_m;
}
