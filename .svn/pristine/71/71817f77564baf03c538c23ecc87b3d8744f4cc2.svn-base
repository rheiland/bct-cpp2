#include <cmath>

#include "bct.h"

/*
 * Preserves a given proportion of the strongest weights in a directed graph.
 * All other weights, as well as those on the main diagonal, are set to zero.
 */
MATRIX_T* BCT_NAMESPACE::threshold_proportional_dir(const MATRIX_T* W, FP_T p) {
	if (safe_mode) check_status(W, SQUARE | DIRECTED, "threshold_proportional_dir");
	
	// n=size(W,1);
	int n = W->size1;
	
	// W(1:n+1:end)=0;
	MATRIX_T* W_thr = copy(W);
	for (int i = 0; i < (int)W_thr->size1; i++) {
		MATRIX_ID(set)(W_thr, i, i, 0.0);
	}
	
	// ind=find(W);
	VECTOR_T* ind = find(W_thr);
	
	// E=sortrows([ind W(ind)], -2);
	VECTOR_T* W_thr_ind = ordinal_index(W_thr, ind);
	VECTOR_T* sort_ind;
	VECTOR_T* sort_W_thr_ind = sort(W_thr_ind, "descend", &sort_ind);
	VECTOR_ID(free)(W_thr_ind);
	VECTOR_ID(free)(sort_W_thr_ind);
	
	// en=round((n^2-n)*p);
	int en = (int)std::floor(n * (n - 1) * p + 0.5);
	
	// W(E(en+1:end,1))=0;
	for (int i = en; i < (int)sort_ind->size; i++) {
		int index = VECTOR_ID(get)(ind, VECTOR_ID(get)(sort_ind, i));
		ordinal_index_assign(W_thr, index, 0.0);
	}
	
	VECTOR_ID(free)(ind);
	VECTOR_ID(free)(sort_ind);
	return W_thr;
}

/*
 * Preserves a given proportion of the strongest weights in an undirected graph.
 * All other weights, as well as those on the main diagonal, are set to zero.
 */
MATRIX_T* BCT_NAMESPACE::threshold_proportional_und(const MATRIX_T* W, FP_T p) {
	if (safe_mode) check_status(W, SQUARE | UNDIRECTED, "threshold_proportional_und");
	
	// n=size(W,1);
	int n = W->size1;
	
	// W(1:n+1:end)=0;
	MATRIX_T* W_thr = triu(W, 1);
	
	// ind=find(W);
	VECTOR_T* ind = find(W_thr);
	
	// E=sortrows([ind W(ind)], -2);
	VECTOR_T* W_thr_ind = ordinal_index(W_thr, ind);
	VECTOR_T* sort_ind;
	VECTOR_T* sort_W_thr_ind = sort(W_thr_ind, "descend", &sort_ind);
	VECTOR_ID(free)(W_thr_ind);
	VECTOR_ID(free)(sort_W_thr_ind);
	
	// en=round((n^2-n)*p);
	int en = (int)std::floor(0.5 * n * (n - 1) * p + 0.5);
	
	// W(E(en+1:end,1))=0;
	for (int i = en; i < (int)sort_ind->size; i++) {
		int index = VECTOR_ID(get)(ind, VECTOR_ID(get)(sort_ind, i));
		ordinal_index_assign(W_thr, index, 0.0);
	}
	for (int i = 0; i < (int)W_thr->size1; i++) {
		for (int j = i + 1; j < (int)W_thr->size2; j++) {
			FP_T value = MATRIX_ID(get)(W_thr, i, j);
			MATRIX_ID(set)(W_thr, j, i, value);
		}
	}
	
	VECTOR_ID(free)(ind);
	VECTOR_ID(free)(sort_ind);
	return W_thr;
}
