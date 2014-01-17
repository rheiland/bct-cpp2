#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Computes the range for each edge (i.e., the shortest path length between the
 * nodes it connects after the edge has been removed from the graph)
 */
MATRIX_T* BCT_NAMESPACE::erange(const MATRIX_T* CIJ, FP_T* eta, MATRIX_T** Eshort, FP_T* fs) {
	if (safe_mode) check_status(CIJ, SQUARE, "erange");
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// K = length(nonzeros(CIJ));
	int k = nnz(CIJ);
	
	// Erange = zeros(N,N);
	MATRIX_T* Erange = zeros(N, N);
	
	// [i,j] = find(CIJ==1);
	MATRIX_T* CIJ_eq_1 = compare_elements(CIJ, fp_equal, 1.0);
	MATRIX_T* find_CIJ_eq_1 = find_ij(CIJ_eq_1);
	MATRIX_ID(free)(CIJ_eq_1);	
	VECTOR_ID(view) i = MATRIX_ID(column)(find_CIJ_eq_1, 0);
	VECTOR_ID(view) j = MATRIX_ID(column)(find_CIJ_eq_1, 1);
	
	// for c=1:length(i)
	for (int c = 0; c < length(&i.vector); c++) {
		
		// CIJcut = CIJ;
		MATRIX_T* CIJcut = copy(CIJ);
		
		// CIJcut(i(c),j(c)) = 0;
		int i_c = (int)VECTOR_ID(get)(&i.vector, c);
		int j_c = (int)VECTOR_ID(get)(&j.vector, c);
		MATRIX_ID(set)(CIJcut, i_c, j_c, 0.0);
		
		// [R,D] = reachdist(CIJcut);
		MATRIX_T* D;
		MATRIX_T* R = reachdist(CIJcut, &D);
		MATRIX_ID(free)(CIJcut);
		MATRIX_ID(free)(R);
		
		// Erange(i(c),j(c)) = D(i(c),j(c))
		MATRIX_ID(set)(Erange, i_c, j_c, MATRIX_ID(get)(D, i_c, j_c));
		MATRIX_ID(free)(D);
	}
	
	MATRIX_ID(free)(find_CIJ_eq_1);
	
	// eta = sum(Erange((Erange>0)&(Erange<Inf)))/length(Erange((Erange>0)&(Erange<Inf)));
	if (eta != NULL) {
		MATRIX_T* Erange_gt_0 = compare_elements(Erange, fp_greater, 0.0);
		MATRIX_T* Erange_lt_inf = compare_elements(Erange, fp_less, GSL_POSINF);
		MATRIX_T* Erange_gt_0_and_Erange_lt_inf = logical_and(Erange_gt_0, Erange_lt_inf);
		VECTOR_T* Erange_idx = logical_index(Erange, Erange_gt_0_and_Erange_lt_inf);
		MATRIX_ID(free)(Erange_gt_0);
		MATRIX_ID(free)(Erange_lt_inf);
		MATRIX_ID(free)(Erange_gt_0_and_Erange_lt_inf);
		*eta = sum(Erange_idx) / (FP_T)length(Erange_idx);
		VECTOR_ID(free)(Erange_idx);
	}
	
	// Eshort = Erange>2;
	MATRIX_T* _Eshort = compare_elements(Erange, fp_greater, 2.0);
	
	// fs = length(nonzeros(Eshort))/K;
	if (fs != NULL) {
		*fs = (FP_T)nnz(_Eshort) / (FP_T)k;
	}
	
	if (Eshort != NULL) *Eshort = _Eshort; else MATRIX_ID(free)(_Eshort);
	return Erange;
}
