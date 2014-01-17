#include "bct.h"

FP_T assortativity(const VECTOR_T*, const MATRIX_T*);

/*
 * Computes assortativity for a directed graph.  Connection weights are ignored.
 */
FP_T BCT_NAMESPACE::assortativity_dir(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE | DIRECTED, "assortativity_dir");
	
	// [id,od,deg] = degrees_dir(CIJ);
	VECTOR_T* deg = degrees_dir(CIJ);
	
	// [i,j] = find(CIJ>0);
	MATRIX_T* CIJ_gt_0 = compare_elements(CIJ, fp_greater, 0.0);
	MATRIX_T* CIJ_gt_0_ij = find_ij(CIJ_gt_0);
	MATRIX_ID(free)(CIJ_gt_0);
	
	FP_T ret = assortativity(deg, CIJ_gt_0_ij);
	VECTOR_ID(free)(deg);
	MATRIX_ID(free)(CIJ_gt_0_ij);
	return ret;
}

/*
 * Computes assortativity for an undirected graph.  Connection weights are
 * ignored.
 */
FP_T BCT_NAMESPACE::assortativity_und(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE | UNDIRECTED, "assortativity_und");
	
	// [deg] = degrees_und(m);
	VECTOR_T* deg = degrees_und(CIJ);
	
	// [i,j] = find(triu(CIJ,1)>0);
	MATRIX_T* triu_CIJ = triu(CIJ, 1);
	MATRIX_T* triu_CIJ_gt_0 = compare_elements(triu_CIJ, fp_greater, 0.0);
	MATRIX_ID(free)(triu_CIJ);
	MATRIX_T* triu_CIJ_gt_0_ij = find_ij(triu_CIJ_gt_0);
	MATRIX_ID(free)(triu_CIJ_gt_0);
	
	FP_T ret = assortativity(deg, triu_CIJ_gt_0_ij);
	VECTOR_ID(free)(deg);
	MATRIX_ID(free)(triu_CIJ_gt_0_ij);
	return ret;
}

FP_T assortativity(const VECTOR_T* deg, const MATRIX_T* ij) {
	using namespace BCT_NAMESPACE;
	
	VECTOR_ID(const_view) i = MATRIX_ID(const_column)(ij, 0);
	VECTOR_ID(const_view) j = MATRIX_ID(const_column)(ij, 1);
	VECTOR_T* degi = VECTOR_ID(alloc)(ij->size1);
	VECTOR_T* degj = VECTOR_ID(alloc)(ij->size1);		
	
	// K = length(i);
	int K = length(&i.vector);
	
	// for k=1:K
	for (int k = 0; k < K; k++) {
		
		// degi(k) = deg(i(k));
		int i_k = (int)VECTOR_ID(get)(&i.vector, k);
		VECTOR_ID(set)(degi, k, VECTOR_ID(get)(deg, i_k));
		
		// degj(k) = deg(j(k));
		int j_k = (int)VECTOR_ID(get)(&j.vector, k);
		VECTOR_ID(set)(degj, k, VECTOR_ID(get)(deg, j_k));
	}
	
	// r = (sum(degi.*degj)/K - (sum(0.5*(degi+degj))/K)^2)/(sum(0.5*(degi.^2+degj.^2))/K - (sum(0.5*(degi+degj))/K)^2);
	
	VECTOR_T* degi_mul_degj = copy(degi);
	VECTOR_ID(mul)(degi_mul_degj, degj);
	FP_T r1 = sum(degi_mul_degj) / (FP_T)K;
	VECTOR_ID(free)(degi_mul_degj);
	
	VECTOR_T* degi_add_degj = copy(degi);
	VECTOR_ID(add)(degi_add_degj, degj);
	FP_T r2 = 0.5 * sum(degi_add_degj) / (FP_T)K;
	VECTOR_ID(free)(degi_add_degj);
	r2 *= r2;
	
	VECTOR_T* degi_pow_2_add_degj_pow_2 = pow_elements(degi, 2);
	VECTOR_T* degj_pow_2 = pow_elements(degj, 2);
	VECTOR_ID(add)(degi_pow_2_add_degj_pow_2, degj_pow_2);
	VECTOR_ID(free)(degj_pow_2);
	FP_T r3 = 0.5 * sum(degi_pow_2_add_degj_pow_2) / (FP_T)K;
	VECTOR_ID(free)(degi_pow_2_add_degj_pow_2);
	
	VECTOR_ID(free)(degi);
	VECTOR_ID(free)(degj);
	
	return (r1 - r2) / (r3 - r2);
}
