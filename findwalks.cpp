#include "bct.h"

/*
 * Finds walks.  Note that there is no twalk argument as its value may overflow
 * a C++ long.  Wq (the main return) and wlq are indexed by path length.  They
 * therefore contain no valid data at index 0.
 */
std::vector<MATRIX_T*> BCT_NAMESPACE::findwalks(const MATRIX_T* CIJ, VECTOR_T** wlq) {
	if (safe_mode) check_status(CIJ, SQUARE, "findwalks");
	
	// CIJ = double(CIJ~=0);
	MATRIX_T* _CIJ = compare_elements(CIJ, fp_not_equal, 0.0);
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// Wq = zeros(N,N,N);
	std::vector<MATRIX_T*> Wq(N + 1);
	Wq[0] = NULL;
	for (int i = 1; i <= N; i++) {
		Wq[i] = zeros(N, N);
	}
	
	// CIJpwr = CIJ;
	MATRIX_T* _CIJpwr = copy(_CIJ);
	
	// Wq(:,:,1) = CIJ;
	Wq[1] = copy(_CIJ);
	
	// for q=2:N
	for (int q = 2; q <= N; q++) {
		
		// CIJpwr = CIJpwr*CIJ;
		MATRIX_T* temp = mul(_CIJpwr, _CIJ);
		MATRIX_ID(free)(_CIJpwr);
		_CIJpwr = temp;
		
		// Wq(:,:,q) = CIJpwr;
		Wq[q] = copy(_CIJpwr);
	}
	
	MATRIX_ID(free)(_CIJ);
	MATRIX_ID(free)(_CIJpwr);
	
	// twalk = sum(sum(sum(Wq)));
	// wlq = reshape(sum(sum(Wq)),1,N);
	if (wlq != NULL) {
		*wlq = VECTOR_ID(alloc)(N + 1);
		VECTOR_ID(set)(*wlq, 0, 0.0);
		for (int i = 1; i <= N; i++) {
			VECTOR_T* sum_Wq_i = sum(Wq[i]);
			VECTOR_ID(set)(*wlq, i, sum(sum_Wq_i));
			VECTOR_ID(free)(sum_Wq_i);
		}
	}
	
	return Wq;
}
