#include <cmath>

#include "bct.h"

/*
 * Generates a random directed binary graph with a hierarchical (fractal)
 * cluster organization.  Cluster connection density starts at 1 and decays as
 * (1 / (E ^ n)), where n is the hierarchical level index.  Cluster size is
 * given by (2 ^ sz_cl).
 */
MATRIX_T* BCT_NAMESPACE::makefractalCIJ(int mx_lvl, FP_T E, int sz_cl, int* K) {
	
	// t = ones(2).*2;
	MATRIX_T* t = MATRIX_ID(alloc)(2, 2);
	MATRIX_ID(set_all)(t, 2.0);
	
	// N = 2^mx_lvl;
	int N = (int)std::pow(2.0, mx_lvl);
	
	// sz_cl = sz_cl-1;
	sz_cl--;
	
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
	
	// ee = mx_lvl-CIJ-sz_cl;
	MATRIX_T* ee = copy(CIJ);
	MATRIX_ID(scale)(ee, -1.0);
	MATRIX_ID(add_constant)(ee, (FP_T)(mx_lvl - sz_cl));
	
	// ee = (ee>0).*ee;
	MATRIX_T* temp = compare_elements(ee, fp_greater, 0.0);
	MATRIX_ID(mul_elements)(temp, ee);
	MATRIX_ID(free)(ee);
	ee = temp;
	
	// prob = (1./(E.^ee)).*(ones(s,s)-eye(s));
	MATRIX_T* E_m = MATRIX_ID(alloc)(s, s);
	MATRIX_ID(set_all)(E_m, E);
	MATRIX_T* neg_ee = copy(ee);
	MATRIX_ID(scale)(neg_ee, -1.0);
	MATRIX_T* prob = pow_elements(E_m, neg_ee);
	MATRIX_ID(free)(E_m);
	MATRIX_ID(free)(neg_ee);
	MATRIX_T* ones_s_sub_eye_s = ones(s, s);
	MATRIX_T* eye_s = eye(s);
	MATRIX_ID(sub)(ones_s_sub_eye_s, eye_s);
	MATRIX_ID(free)(eye_s);
	MATRIX_ID(mul_elements)(prob, ones_s_sub_eye_s);
	MATRIX_ID(free)(ones_s_sub_eye_s);
	
	// CIJ = (prob>rand(N));
	MATRIX_T* rand_N = rand(N);
	MATRIX_ID(free)(CIJ);
	CIJ = compare_elements(prob, fp_greater, rand_N);
	MATRIX_ID(free)(rand_N);
	
	// K = sum(sum(CIJ));
	if (K != NULL) {
		VECTOR_T* sum_CIJ = sum(CIJ);
		*K = sum(sum_CIJ);
		VECTOR_ID(free)(sum_CIJ);
	}
	
	return CIJ;
}
