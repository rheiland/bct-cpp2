#include <cmath>

#include "bct.h"

/*
 * Generates a random directed binary graph with equal-sized clusters placed on
 * the diagonal and the remaining connections distributed randomly among them.
 * N must be a power of 2, and the cluster size is given by (2 ^ sz_cl).
 */
MATRIX_T* BCT_NAMESPACE::makeevenCIJ(int N, int K, int sz_cl) {

	// mx_lvl = floor(log2(N));
	int mx_lvl = (int)std::floor(std::log((FP_T)N) / std::log(2.0));
	
	// sz_cl = sz_cl-1;
	sz_cl--;

	// t = ones(2).*2;
	MATRIX_T* t = MATRIX_ID(alloc)(2, 2);
	MATRIX_ID(set_all)(t, 2.0);

	// Nlvl = 2^mx_lvl;
	FP_T Nlvl = std::pow(2.0, mx_lvl);
	
	// N = Nlvl;
	N = (int)Nlvl;
	
	MATRIX_T* CIJ = NULL;
	
	// for lvl=1:mx_lvl-1
	for (int lvl = 1; lvl <= mx_lvl - 1; lvl++) {
		
		// CIJ = ones(2^(lvl+1),2^(lvl+1));
		int n = (int)std::pow(2.0, lvl + 1);
		CIJ = ones(n, n);
		
		// group1 = [1:size(CIJ,1)/2];
		VECTOR_T* group1 = sequence(0, n / 2 - 1);
		
		// group2 = [size(CIJ,1)/2+1:size(CIJ,1)];
		VECTOR_T* group2 = sequence(n / 2, n - 1);
		
		// CIJ(group1,group1) = t;
		ordinal_index_assign(CIJ, group1, group1, t);
		VECTOR_ID(free)(group1);
		
		// CIJ(group2,group2) = t;
		ordinal_index_assign(CIJ, group2, group2, t);
		VECTOR_ID(free)(group2);
		
		// CIJ = CIJ+ones(size(CIJ,1),size(CIJ,1));
		MATRIX_ID(add_constant)(CIJ, 1.0);
		
		// t = CIJ;
		MATRIX_ID(free)(t);
		t = CIJ;
	}
	
	// s = size(CIJ,1);
	int s = CIJ->size1;
	
	// CIJ = CIJ-ones(s,s)-mx_lvl.*eye(s);
	MATRIX_ID(add_constant)(CIJ, -1.0);
	MATRIX_T* mx_lvl_mul_eye_s = eye(s);
	MATRIX_ID(scale)(mx_lvl_mul_eye_s, (FP_T)mx_lvl);
	MATRIX_ID(sub)(CIJ, mx_lvl_mul_eye_s);
	MATRIX_ID(free)(mx_lvl_mul_eye_s);

	// CIJp = (CIJ>=(mx_lvl-sz_cl));
	MATRIX_T* CIJp = compare_elements(CIJ, fp_greater_or_equal, (FP_T)(mx_lvl - sz_cl));
	MATRIX_ID(free)(CIJ);

	// CIJc = (CIJp==1);
	MATRIX_T* CIJc = compare_elements(CIJp, fp_equal, 1.0);
	MATRIX_ID(free)(CIJp);
	
	// remK = K-nnz(CIJc);
	int remK = K - nnz(CIJc);
	
	if (remK > 0) {
		
		// [a,b] = find(~(CIJc+eye(N)));
		MATRIX_T* CIJc_add_eye_N = copy(CIJc);
		MATRIX_T* eye_N = eye(N);
		MATRIX_ID(add)(CIJc_add_eye_N, eye_N);
		MATRIX_ID(free)(eye_N);
		MATRIX_T* not_CIJc_add_eye_N = logical_not(CIJc_add_eye_N);
		MATRIX_ID(free)(CIJc_add_eye_N);
		MATRIX_T* find_not_CIJc_add_eye_N = find_ij(not_CIJc_add_eye_N);
		MATRIX_ID(free)(not_CIJc_add_eye_N);
		VECTOR_ID(view) a = MATRIX_ID(column)(find_not_CIJc_add_eye_N, 0);
		VECTOR_ID(view) b = MATRIX_ID(column)(find_not_CIJc_add_eye_N, 1);
		
		// rp = randperm(length(a));
		gsl_permutation* rp = randperm(length(&a.vector));
		VECTOR_T* rp_v = to_vector(rp);
		gsl_permutation_free(rp);
		VECTOR_ID(view) rp_subv = VECTOR_ID(subvector)(rp_v, 0, remK);
		
		// a = a(rp(1:remK));
		VECTOR_T* a_rp_subv = ordinal_index(&a.vector, &rp_subv.vector);
		VECTOR_T* ab_indices = sequence(0, remK - 1);
		ordinal_index_assign(&a.vector, ab_indices, a_rp_subv);
		VECTOR_ID(free)(a_rp_subv);
		
		// b = b(rp(1:remK));
		VECTOR_T* b_rp_subv = ordinal_index(&b.vector, &rp_subv.vector);
		VECTOR_ID(free)(rp_v);
		ordinal_index_assign(&b.vector, ab_indices, b_rp_subv);
		VECTOR_ID(free)(ab_indices);
		VECTOR_ID(free)(b_rp_subv);
		
		// for i=1:remK
		for (int i = 0; i < remK; i++) {
			
			// CIJc(a(i),b(i)) = 1;
			MATRIX_ID(set)(CIJc, (int)VECTOR_ID(get)(&a.vector, i), (int)VECTOR_ID(get)(&b.vector, i), 1.0);
		}
		
		MATRIX_ID(free)(find_not_CIJc_add_eye_N);
	}

	// CIJ = CIJc;
	return CIJc;
}
