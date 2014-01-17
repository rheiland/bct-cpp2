#include "bct.h"

FP_T matching_ind(const VECTOR_T*, const VECTOR_T*, int, int, int);

/*
 * Computes matching index for all connections.
 */
MATRIX_T* BCT_NAMESPACE::matching_ind(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE, "matching_ind");
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// Mall = zeros(N,N);
	MATRIX_T* Mall = zeros(N, N);
	
	// for i=1:N-1
	for (int i = 0; i < N - 1; i++) {
		
		// c1 = [CIJ(:,i)' CIJ(i,:)];
		VECTOR_ID(const_view) CIJ_col_i = MATRIX_ID(const_column)(CIJ, i);
		VECTOR_ID(const_view) CIJ_row_i = MATRIX_ID(const_row)(CIJ, i);
		VECTOR_T* c1 = concatenate(&CIJ_col_i.vector, &CIJ_row_i.vector);
		
		// for j=i+1:N
		for (int j = i + 1; j < N; j++) {
			
			// c2 = [CIJ(:,j)' CIJ(j,:)];
			VECTOR_ID(const_view) CIJ_col_j = MATRIX_ID(const_column)(CIJ, j);
			VECTOR_ID(const_view) CIJ_row_j = MATRIX_ID(const_row)(CIJ, j);
			VECTOR_T* c2 = concatenate(&CIJ_col_j.vector, &CIJ_row_j.vector);
			
			MATRIX_ID(set)(Mall, i, j, matching_ind(c1, c2, i, j, N));
			VECTOR_ID(free)(c2);
		}
		
		VECTOR_ID(free)(c1);
	}
	
	return Mall;
}

/*
 * Computes matching index for incoming connections.
 */
MATRIX_T* BCT_NAMESPACE::matching_ind_in(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE, "matching_ind_in");
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// Min = zeros(N,N);
	MATRIX_T* Min = zeros(N, N);
	
	// for i=1:N-1
	for (int i = 0; i < N - 1; i++) {
		
		// c1 = CIJ(:,i);
		VECTOR_ID(const_view) c1 = MATRIX_ID(const_column)(CIJ, i);
		
		// for j=i+1:N
		for (int j = i + 1; j < N; j++) {
			
			// c2 = CIJ(:,j);
			VECTOR_ID(const_view) c2 = MATRIX_ID(const_column)(CIJ, j);
			
			MATRIX_ID(set)(Min, i, j, matching_ind(&c1.vector, &c2.vector, i, j, 0));
		}
	}
	
	return Min;
}

/*
 * Computes matching index for outgoing connections.
 */
MATRIX_T* BCT_NAMESPACE::matching_ind_out(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE, "matching_ind_out");
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// Mout = zeros(N,N);
	MATRIX_T* Mout = zeros(N, N);
	
	// for i=1:N-1
	for (int i = 0; i < N - 1; i++) {
		
		// c1 = CIJ(:,i);
		VECTOR_ID(const_view) c1 = MATRIX_ID(const_row)(CIJ, i);
		
		// for j=i+1:N
		for (int j = i + 1; j < N; j++) {
			
			// c2 = CIJ(:,j);
			VECTOR_ID(const_view) c2 = MATRIX_ID(const_row)(CIJ, j);
			
			MATRIX_ID(set)(Mout, i, j, matching_ind(&c1.vector, &c2.vector, i, j, 0));
		}
	}
	
	return Mout;
}

FP_T matching_ind(const VECTOR_T* c1, const VECTOR_T* c2, int i, int j, int N) {
	using namespace BCT_NAMESPACE;
	
	// use = ~(~c1&~c2);
	VECTOR_T* use = logical_or(c1, c2);
	
	// use(i) = 0;  use(i+N) = 0;
	VECTOR_ID(set)(use, i, 0.0);
	VECTOR_ID(set)(use, i + N, 0.0);
	
	// use(j) = 0;  use(j+N) = 0;
	VECTOR_ID(set)(use, j, 0.0);
	VECTOR_ID(set)(use, j + N, 0.0);
	
	// ncon = sum(c1(use))+sum(c2(use));
	VECTOR_T* c1_use = logical_index(c1, use);
	VECTOR_T* c2_use = logical_index(c2, use);
	VECTOR_ID(free)(use);
	FP_T ncon = 0.0;
	if (c1_use != NULL) {
		ncon += sum(c1_use);
	}
	if (c2_use != NULL) {
		ncon += sum(c2_use);
	}
	
	// if (ncon==0)
	FP_T ret;
	if (fp_zero(ncon)) {
		
		// Mall(i,j) = 0;
		ret = 0.0;
	} else {
		
		// Mall(i,j) = 2*(sum(c1(use)&c2(use))/ncon);
		if (c1_use == NULL || c2_use == NULL) {
			ret = 0.0;
		} else {
			VECTOR_T* log_and = logical_and(c1_use, c2_use);
			ret = 2.0 * sum(log_and) / ncon;
			VECTOR_ID(free)(log_and);
		}
	}
	
	if (c1_use != NULL) VECTOR_ID(free)(c1_use);
	if (c2_use != NULL) VECTOR_ID(free)(c2_use);
	return ret;
}
