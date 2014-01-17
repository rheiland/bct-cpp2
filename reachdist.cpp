#include <gsl/gsl_math.h>

#include "bct.h"

void reachdist2(const MATRIX_T*, MATRIX_T**, MATRIX_T**, MATRIX_T**, int, int*, VECTOR_T*, VECTOR_T*);

/*
 * Computes reachability and distance matrices based on the power of the
 * adjacency matrix.
 */
MATRIX_T* BCT_NAMESPACE::reachdist(const MATRIX_T* CIJ, MATRIX_T** D) {
	if (safe_mode) check_status(CIJ, SQUARE, "reachdist");
	
	// R = CIJ;
	MATRIX_T* R = copy(CIJ);
	
	// D = CIJ;
	if (D != NULL) {
		*D = copy(CIJ);
	}
	
	// powr = 2;
	int powr = 2;
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// CIJpwr = CIJ;
	MATRIX_T* CIJpwr = copy(CIJ);
	
	// id = sum(CIJ,1);
	VECTOR_T* id = sum(CIJ, 1);
	
	// od = sum(CIJ,2)';
	VECTOR_T* od = sum(CIJ, 2);
	
	// id_0 = find(id==0);
	VECTOR_T* id_eq_0 = compare_elements(id, fp_equal, 0.0);
	VECTOR_ID(free)(id);
	VECTOR_T* id_0 = find(id_eq_0);
	VECTOR_ID(free)(id_eq_0);
	
	// od_0 = find(od==0);
	VECTOR_T* od_eq_0 = compare_elements(od, fp_equal, 0.0);
	VECTOR_ID(free)(od);
	VECTOR_T* od_0 = find(od_eq_0);
	VECTOR_ID(free)(od_eq_0);
	
	// col = setxor(1:N,id_0);
	VECTOR_T* i_all = sequence(0, N - 1);
	VECTOR_T* col = setxor(i_all, id_0);
	
	// row = setxor(1:N,od_0);
	VECTOR_T* row = setxor(i_all, od_0);
	VECTOR_ID(free)(i_all);
	
	// [R,D,powr] = reachdist2(CIJ,CIJpwr,R,D,N,powr,col,row);
	reachdist2(CIJ, &CIJpwr, &R, D, N, &powr, col, row);
	MATRIX_ID(free)(CIJpwr);
	if (col != NULL) {
		VECTOR_ID(free)(col);
	}
	if (row != NULL) {
		VECTOR_ID(free)(row);
	}
	
	if (D != NULL) {
		
		// D = powr - D+1;
		MATRIX_ID(scale)(*D, -1.0);
		MATRIX_ID(add_constant)(*D, (FP_T)(powr + 1));
		
		// D(D==(N+2)) = Inf;
		MATRIX_T* D_eq_N_add_2 = compare_elements(*D, fp_equal, (FP_T)(N + 2));
		logical_index_assign(*D, D_eq_N_add_2, GSL_POSINF);
		MATRIX_ID(free)(D_eq_N_add_2);
		
		// D(:,id_0) = Inf;
		VECTOR_T* D_rows_cols = sequence(0, N - 1);
		if (id_0 != NULL) {
			ordinal_index_assign(*D, D_rows_cols, id_0, GSL_POSINF);
		}
		
		// D(od_0,:) = Inf;
		if (od_0 != NULL) {
			ordinal_index_assign(*D, od_0, D_rows_cols, GSL_POSINF);
		}
		VECTOR_ID(free)(D_rows_cols);
	}
	
	if (id_0 != NULL) {
		VECTOR_ID(free)(id_0);
	}
	if (od_0 != NULL) {
		VECTOR_ID(free)(od_0);
	}
	return R;
}

void reachdist2(const MATRIX_T* CIJ, MATRIX_T** CIJpwr, MATRIX_T** R, MATRIX_T** D, int N, int* powr, VECTOR_T* col, VECTOR_T* row) {
	using namespace BCT_NAMESPACE;
	
	// CIJpwr = CIJpwr*CIJ;
	MATRIX_T* temp = mul(*CIJpwr, CIJ);
	MATRIX_ID(free)(*CIJpwr);
	*CIJpwr = temp;
	
	// R = double(R | ((CIJpwr)~=0));
	MATRIX_T* CIJpwr_neq_0 = compare_elements(*CIJpwr, fp_not_equal, 0.0);
	MATRIX_T* R_or_CIJpwr_neq_0 = logical_or(*R, CIJpwr_neq_0);
	MATRIX_ID(free)(CIJpwr_neq_0);
	MATRIX_ID(free)(*R);
	*R = R_or_CIJpwr_neq_0;
	
	// D = D+R;
	if (D != NULL) {
		MATRIX_ID(add)(*D, *R);
	}
	
	// if ((powr<=N)&&(~isempty(nonzeros(R(row,col)==0))))
	if (*powr <= N) {
		if (row != NULL && col != NULL) {
			MATRIX_T* R_idx = ordinal_index(*R, row, col);
			MATRIX_T* R_idx_eq_0 = compare_elements(R_idx, fp_equal, 0.0);
			MATRIX_ID(free)(R_idx);
			VECTOR_T* nonzeros_R_idx_eq_0 = nonzeros(R_idx_eq_0);
			MATRIX_ID(free)(R_idx_eq_0);
			if (nonzeros_R_idx_eq_0 != NULL) {
				VECTOR_ID(free)(nonzeros_R_idx_eq_0);
				
				// powr = powr+1;
				(*powr)++;
				
				// [R,D,powr] = reachdist2(CIJ,CIJpwr,R,D,N,powr,col,row);
				reachdist2(CIJ, CIJpwr, R, D, N, powr, col, row);
			}
		}
	}
}
