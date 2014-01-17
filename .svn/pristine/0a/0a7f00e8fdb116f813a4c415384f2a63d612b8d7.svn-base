#include "bct.h"

/*
 * Computes the joint degree distribution, a matrix in which the value of each
 * element (u, v) is the number of nodes with u outgoing connections and v
 * incoming connections.
 */
MATRIX_T* BCT_NAMESPACE::jdegree(const MATRIX_T* CIJ) {
	if (safe_mode) check_status(CIJ, SQUARE, "jdegree");
	
	// CIJ = double(CIJ~=0);
	MATRIX_T* _CIJ = compare_elements(CIJ, fp_not_equal, 0.0);
	
	// N = size(CIJ,1);
	int N = _CIJ->size1;
	
	// id = sum(CIJ,1);
	VECTOR_T* id = sum(_CIJ, 1);
	
	// od = sum(CIJ,2)';
	VECTOR_T* od = sum(_CIJ, 2);
	MATRIX_ID(free)(_CIJ);
	
	// szJ = max(max(id,od))+1;
	FP_T max_id = VECTOR_ID(max)(id);
	FP_T max_od = VECTOR_ID(max)(od);
	int szJ = (int)max(max_id, max_od) + 1;
	
	// J = zeros(szJ);
	MATRIX_T* J = zeros(szJ);
	
	// for i=1:N
	for (int i = 0; i < N; i++) {
		
		// J(id(i)+1,od(i)+1) = J(id(i)+1,od(i)+1) + 1;
		int id_i = (int)VECTOR_ID(get)(id, i);
		int od_i = (int)VECTOR_ID(get)(od, i);
		MATRIX_ID(set)(J, id_i, od_i, MATRIX_ID(get)(J, id_i, od_i) + 1.0);
	}
	
	VECTOR_ID(free)(id);
	VECTOR_ID(free)(od);
	return J;
}

/*
 * Given a joint degree distribution matrix, returns the number of nodes with
 * in-degree = out-degree.
 */
int BCT_NAMESPACE::jdegree_bl(const MATRIX_T* J) {
	
	// J_bl = sum(diag(J));
	VECTOR_ID(const_view) diag_J = MATRIX_ID(const_diagonal)(J);
	return (int)sum(&diag_J.vector);
}

/*
 * Given a joint degree distribution matrix, returns the number of nodes with
 * in-degree > out-degree.
 */
int BCT_NAMESPACE::jdegree_id(const MATRIX_T* J) {
	
	// J_id = sum(sum(tril(J,-1)));
	MATRIX_T* tril_J = tril(J, -1);
	VECTOR_T* sum_tril_J = sum(tril_J);
	int J_id = (int)sum(sum_tril_J);
	MATRIX_ID(free)(tril_J);
	VECTOR_ID(free)(sum_tril_J);
	return J_id;
}

/*
 * Given a joint degree distribution matrix, returns the number of nodes with
 * out-degree > in-degree.
 */
int BCT_NAMESPACE::jdegree_od(const MATRIX_T* J) {
	
	// J_od = sum(sum(triu(J,1)));
	MATRIX_T* triu_J = triu(J, 1);
	VECTOR_T* sum_triu_J = sum(triu_J);
	int J_od = (int)sum(sum_triu_J);
	MATRIX_ID(free)(triu_J);
	VECTOR_ID(free)(sum_triu_J);	
	return J_od;
}
