#include "bct.h"

/*
 * Generates a random directed binary graph with a Toeplitz organization.
 */
MATRIX_T* BCT_NAMESPACE::maketoeplitzCIJ(int N, int K, FP_T s) {
	
	// profile = normpdf([1:N-1],0.5,s);
	VECTOR_T* indices = sequence(1, N - 1);
	VECTOR_T* profile = normpdf(indices, 0.5, s);
	VECTOR_ID(free)(indices);
	
	// template = toeplitz([0 profile],[0 profile]);
	VECTOR_T* temp = concatenate(0.0, profile);
	VECTOR_ID(free)(profile);
	profile = temp;
	MATRIX_T* _template = toeplitz(profile, profile);
	VECTOR_ID(free)(profile);
	
	// template = template.*(K./sum(sum(template)))
	VECTOR_T* sum__template = sum(_template);
	FP_T sum_sum__template = sum(sum__template);
	VECTOR_ID(free)(sum__template);
	MATRIX_ID(scale)(_template, (FP_T)K / sum_sum__template);
	
	// CIJ = zeros(N);
	MATRIX_T* CIJ = zeros(N);
	
	// while ((sum(sum(CIJ)) ~= K))
	VECTOR_T* sum_CIJ = sum(CIJ);
	FP_T sum_sum_CIJ = sum(sum_CIJ);
	VECTOR_ID(free)(sum_CIJ);
	while ((int)sum_sum_CIJ != K) {
		
		// CIJ = (rand(N)<template);
		MATRIX_T* rand_N = rand(N);
		MATRIX_ID(free)(CIJ);
		CIJ = compare_elements(rand_N, fp_less, _template);
		MATRIX_ID(free)(rand_N);
		
		sum_CIJ = sum(CIJ);
		sum_sum_CIJ = sum(sum_CIJ);
		VECTOR_ID(free)(sum_CIJ);
	}
	
	MATRIX_ID(free)(_template);
	return CIJ;
}
