#include <cmath>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>

#include "bct.h"

VECTOR_T* modularity(const MATRIX_T* B, int N, FP_T m);

/*
 * Detects communities in a directed graph via Newman modularity.  Since GSL
 * solves eigensystems differently from MATLAB, communities may be numbered
 * differently.
 */
FP_T BCT_NAMESPACE::modularity_dir(const MATRIX_T* A, VECTOR_T** Ci) {
	if (safe_mode) check_status(A, SQUARE | DIRECTED, "modularity_dir");
	
	// N=length(A);
	int N = length(A);
	
	// n_perm = randperm(N);
	gsl_permutation* n_perm = randperm(N);
	
	// A = A(n_perm,n_perm);
	VECTOR_T* n_perm_v = to_vector(n_perm);
	gsl_permutation_free(n_perm);
	MATRIX_T* A_perm = MATRIX_ID(alloc)(N, N);
	ordinal_index_assign(A_perm, n_perm_v, n_perm_v, A);
	
	// Ki=sum(A,1);
	VECTOR_T* Ki = sum(A_perm, 1);
	
	// Ko=sum(A,2);
	VECTOR_T* Ko = sum(A_perm, 2);
	
	// m=sum(Ki);
	FP_T m = sum(Ki);
	
	// b=A-(Ko*Ki).'/m;
	MATRIX_T* Ko_m = to_column_matrix(Ko);
	VECTOR_ID(free)(Ko);
	MATRIX_T* Ki_m = to_row_matrix(Ki);
	VECTOR_ID(free)(Ki);
	MATRIX_T* Ko_mul_Ki_transpose = mul(Ko_m, Ki_m);
	MATRIX_ID(free)(Ko_m);
	MATRIX_ID(free)(Ki_m);
	MATRIX_ID(transpose)(Ko_mul_Ki_transpose);
	MATRIX_ID(scale)(Ko_mul_Ki_transpose, 1.0 / m);
	MATRIX_T* b = copy(A_perm);
	MATRIX_ID(free)(A_perm);
	MATRIX_ID(sub)(b, Ko_mul_Ki_transpose);
	MATRIX_ID(free)(Ko_mul_Ki_transpose);
	
	// B=b+b.';
	MATRIX_T* b_transpose = MATRIX_ID(alloc)(b->size2, b->size1);
	MATRIX_ID(transpose_memcpy)(b_transpose, b);
	MATRIX_T* B = b;
	MATRIX_ID(add)(B, b_transpose);
	MATRIX_ID(free)(b_transpose);
	
	VECTOR_T* _Ci = modularity(B, N, m);
	
	// s=Ci(:,ones(1,N));
	MATRIX_T* s = MATRIX_ID(alloc)(N, N);
	for (int i = 0; i < N; i++) {
		MATRIX_ID(set_col)(s, i, _Ci);
	}
	
	// Q=~(s-s.').*B/(2*m);
	MATRIX_T* s_transpose = MATRIX_ID(alloc)(s->size2, s->size1);
	MATRIX_ID(transpose_memcpy)(s_transpose, s);
	MATRIX_ID(sub)(s, s_transpose);
	MATRIX_ID(free)(s_transpose);
	MATRIX_T* Q_m = logical_not(s);
	MATRIX_ID(free)(s);
	MATRIX_ID(mul_elements)(Q_m, B);
	MATRIX_ID(free)(B);
	MATRIX_ID(scale)(Q_m, 1.0 / (2.0 * m));
	
	// Q=sum(Q(:));
	VECTOR_T* sum_Q_m = sum(Q_m);
	MATRIX_ID(free)(Q_m);
	FP_T Q = sum(sum_Q_m);
	VECTOR_ID(free)(sum_Q_m);
	
	// Ci_corrected = zeros(N,1);
	VECTOR_T* _Ci_corrected = VECTOR_ID(alloc)(N);
	
	// Ci_corrected(n_perm) = Ci;
	for (int i = 0; i < N; i++) {
		int index = (int)VECTOR_ID(get)(n_perm_v, i);
		FP_T value = VECTOR_ID(get)(_Ci, i);
		VECTOR_ID(set)(_Ci_corrected, index, value);
	}
	VECTOR_ID(free)(n_perm_v);
	
	// Ci = Ci_corrected;
	VECTOR_ID(free)(_Ci);
	_Ci = _Ci_corrected;
	
	if (Ci != NULL) *Ci = _Ci; else VECTOR_ID(free)(_Ci);
	return Q;
}

/*
 * Detects communities in an undirected graph via Newman modularity.  Since GSL
 * solves eigensystems differently from MATLAB, communities may be numbered
 * differently.
 */
FP_T BCT_NAMESPACE::modularity_und(const MATRIX_T* A, VECTOR_T** Ci) {
	if (safe_mode) check_status(A, SQUARE | UNDIRECTED, "modularity_und");
	
	// N=length(A);
	int N = length(A);
	
	// n_perm = randperm(N);
	gsl_permutation* n_perm = randperm(N);
	
	// A = A(n_perm,n_perm);
	VECTOR_T* n_perm_v = to_vector(n_perm);
	gsl_permutation_free(n_perm);
	MATRIX_T* A_perm = MATRIX_ID(alloc)(N, N);
	ordinal_index_assign(A_perm, n_perm_v, n_perm_v, A);
	
	// K=sum(A);
	VECTOR_T* K = sum(A_perm);
	
	// m=sum(K);
	FP_T m = sum(K);
	
	// B=A-(K.'*K)/m;
	MATRIX_T* K_m_transpose = to_column_matrix(K);
	MATRIX_T* K_m = to_row_matrix(K);
	VECTOR_ID(free)(K);
	MATRIX_T* K_m_transpose_mul_K_m = mul(K_m_transpose, K_m);
	MATRIX_ID(free)(K_m_transpose);
	MATRIX_ID(free)(K_m);
	MATRIX_ID(scale)(K_m_transpose_mul_K_m, 1.0 / m);
	MATRIX_T* B = copy(A_perm);
	MATRIX_ID(free)(A_perm);
	MATRIX_ID(sub)(B, K_m_transpose_mul_K_m);
	MATRIX_ID(free)(K_m_transpose_mul_K_m);
	
	VECTOR_T* _Ci = modularity(B, N, m);
	
	// s=Ci(:,ones(1,N));
	MATRIX_T* s = MATRIX_ID(alloc)(N, N);
	for (int i = 0; i < N; i++) {
		MATRIX_ID(set_col)(s, i, _Ci);
	}
	
	// Q=~(s-s.').*B/m;
	MATRIX_T* s_transpose = MATRIX_ID(alloc)(s->size2, s->size1);
	MATRIX_ID(transpose_memcpy)(s_transpose, s);
	MATRIX_ID(sub)(s, s_transpose);
	MATRIX_ID(free)(s_transpose);
	MATRIX_T* Q_m = logical_not(s);
	MATRIX_ID(free)(s);
	MATRIX_ID(mul_elements)(Q_m, B);
	MATRIX_ID(free)(B);
	MATRIX_ID(scale)(Q_m, 1.0 / m);
	
	// Q=sum(Q(:));
	VECTOR_T* sum_Q_m = sum(Q_m);
	MATRIX_ID(free)(Q_m);
	FP_T Q = sum(sum_Q_m);
	VECTOR_ID(free)(sum_Q_m);
	
	// Ci_corrected = zeros(N,1);
	VECTOR_T* _Ci_corrected = VECTOR_ID(alloc)(N);
	
	// Ci_corrected(n_perm) = Ci;
	for (int i = 0; i < N; i++) {
		int index = (int)VECTOR_ID(get)(n_perm_v, i);
		FP_T value = VECTOR_ID(get)(_Ci, i);
		VECTOR_ID(set)(_Ci_corrected, index, value);
	}
	VECTOR_ID(free)(n_perm_v);
	
	// Ci = Ci_corrected;
	VECTOR_ID(free)(_Ci);
	_Ci = _Ci_corrected;
	
	if (Ci != NULL) *Ci = _Ci; else VECTOR_ID(free)(_Ci);
	return Q;
}

VECTOR_T* modularity(const MATRIX_T* B, int N, FP_T m) {
	using namespace BCT_NAMESPACE;
	
	// Ci=ones(N,1);
	VECTOR_T* Ci = ones_vector(N);
	
	// cn=1;
	int cn = 1;
	
	// U=[1 0];
	VECTOR_T* U = VECTOR_ID(alloc)(2);
	VECTOR_ID(set)(U, 0, 1.0);
	VECTOR_ID(set)(U, 1, 0.0);
	
	// ind=1:N;
	VECTOR_T* ind = sequence(0, N - 1);
	
	// Bg=B;
	MATRIX_T* Bg = copy(B);
	
	// Ng=N;
	int Ng = N;
	
	gsl_eigen_symmv_workspace* eig = gsl_eigen_symmv_alloc(Bg->size1);
	
	// while U(1)
	while (fp_nonzero(VECTOR_ID(get)(U, 0))) {
		
		// [V D]=eig(Bg);
		// [d1 i1]=max(diag(D));
		// v1=V(:,i1);
		if (eig->size != Bg->size1) {
			gsl_eigen_symmv_free(eig);
			eig = gsl_eigen_symmv_alloc(Bg->size1);
		}
		gsl_matrix* temp = to_matrix_double(Bg);
		gsl_vector* eval = gsl_vector_alloc(Bg->size1);
		gsl_matrix* evec = gsl_matrix_alloc(Bg->size1, Bg->size1);
		gsl_eigen_symmv(temp, eval, evec, eig);
		gsl_matrix_free(temp);
		gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);
		gsl_vector_free(eval);
		gsl_vector* v1_double = gsl_vector_alloc(evec->size1);
		gsl_matrix_get_col(v1_double, evec, 0);
		gsl_matrix_free(evec);
		VECTOR_T* v1 = to_vector(v1_double);
		gsl_vector_free(v1_double);
		
		// S=ones(Ng,1);
		MATRIX_T* S = ones(Ng, 1);
		
		// S(v1<0)=-1;
		VECTOR_T* v1_lt_0 = compare_elements(v1, fp_less, 0.0);
		VECTOR_ID(free)(v1);
		logical_index_assign(S, v1_lt_0, -1.0);
		VECTOR_ID(free)(v1_lt_0);
		
 		// q=S.'*Bg*S;
		MATRIX_T* S_transpose = MATRIX_ID(alloc)(S->size2, S->size1);
		MATRIX_ID(transpose_memcpy)(S_transpose, S);
		MATRIX_T* S_transpose_mul_Bg = mul(S_transpose, Bg);
		MATRIX_ID(free)(S_transpose);
		MATRIX_T* q_m = mul(S_transpose_mul_Bg, S);
		MATRIX_ID(free)(S_transpose_mul_Bg);
		FP_T q = MATRIX_ID(get)(q_m, 0, 0);
		MATRIX_ID(free)(q_m);
		
		// if q>1e-10
		if (q > 1e-10) {
			
			// qmax=q;
			FP_T qmax = q;
			
			// Bg(logical(eye(Ng)))=0;
			MATRIX_T* eye_Ng = eye(Ng);
			logical_index_assign(Bg, eye_Ng, 0.0);
			MATRIX_ID(free)(eye_Ng);
			
			// indg=ones(Ng,1);
			MATRIX_T* indg = ones(Ng, 1);
			
			// Sit=S;
			MATRIX_T* Sit = copy(S);
			
			// while any(indg);
			VECTOR_T* any_indg = any(indg);
			bool any_indg_bool = to_bool(any_indg);
			VECTOR_ID(free)(any_indg);
			while (any_indg_bool) {
				
				// Qit=qmax-4*Sit.*(Bg*Sit);
				MATRIX_T* Qit = mul(Bg, Sit);
				MATRIX_ID(mul_elements)(Qit, Sit);
				MATRIX_ID(scale)(Qit, -4.0);
				MATRIX_ID(add_constant)(Qit, qmax);
				
				// qmax=max(Qit.*indg);
				MATRIX_T* Qit_mul_indg = copy(Qit);
				MATRIX_ID(mul_elements)(Qit_mul_indg, indg);
				VECTOR_T* qmax_v = max(Qit_mul_indg);
				MATRIX_ID(free)(Qit_mul_indg);
				qmax = VECTOR_ID(get)(qmax_v, 0);
				VECTOR_ID(free)(qmax_v);
				
				// TODO: Fix precision issue (differences of 1e-14 for macaque47/71)
				// imax=(Qit==qmax);
				MATRIX_T* imax = compare_elements(Qit, fp_equal, qmax);
				MATRIX_ID(free)(Qit);
				
				// Sit(imax)=-Sit(imax);
				VECTOR_T* neg_Sit_imax = logical_index(Sit, imax);
				if (neg_Sit_imax != NULL) {
					VECTOR_ID(scale)(neg_Sit_imax, -1.0);
					logical_index_assign(Sit, imax, neg_Sit_imax);
					VECTOR_ID(free)(neg_Sit_imax);
				}
				
				// indg(imax)=nan;
				logical_index_assign(indg, imax, GSL_NAN);
				MATRIX_ID(free)(imax);
				
				// if qmax>q;
				if (qmax > q) {
					
					// q=qmax;
					q = qmax;
					
					// S=Sit;
					MATRIX_ID(free)(S);
					S = copy(Sit);
				}
				
				any_indg = any(indg);
				any_indg_bool = to_bool(any_indg);
				VECTOR_ID(free)(any_indg);
			}
			
			// if(abs(sum(S))==Ng
			VECTOR_T* sum_S_v = sum(S);
			FP_T sum_S = sum(sum_S_v);
			VECTOR_ID(free)(sum_S_v);
			if (fp_equal(std::abs(sum_S), (FP_T)Ng)) {
				
				// U(1)=[];
				VECTOR_T* temp = VECTOR_ID(alloc)(U->size - 1);
				VECTOR_ID(view) U_subv = VECTOR_ID(subvector)(U, 1, U->size - 1);
				VECTOR_ID(memcpy)(temp, &U_subv.vector);
				VECTOR_ID(free)(U);
				U = temp;
			} else {
				
				// cn=cn+1;
				cn++;
				
				if (ind != NULL) {
					
					// Ci(ind(S==1))=U(1);
					MATRIX_T* S_eq_1_m = compare_elements(S, fp_equal, 1.0);
					VECTOR_T* S_eq_1 = to_vector(S_eq_1_m);
					MATRIX_ID(free)(S_eq_1_m);
					VECTOR_T* ind_idx = logical_index(ind, S_eq_1);
					VECTOR_ID(free)(S_eq_1);
					if (ind_idx != NULL) {
						ordinal_index_assign(Ci, ind_idx, VECTOR_ID(get)(U, 0));
						VECTOR_ID(free)(ind_idx);
					}
					
					// Ci(ind(S==-1))=cn;
					MATRIX_T* S_eq_neg_1_m = compare_elements(S, fp_equal, -1.0);
					VECTOR_T* S_eq_neg_1 = to_vector(S_eq_neg_1_m);
					MATRIX_ID(free)(S_eq_neg_1_m);
					ind_idx = logical_index(ind, S_eq_neg_1);
					VECTOR_ID(free)(S_eq_neg_1);
					if (ind_idx != NULL) {
						ordinal_index_assign(Ci, ind_idx, (FP_T)cn);
						VECTOR_ID(free)(ind_idx);
					}
				}
				
				// U=[cn U];
				VECTOR_T* temp = concatenate((FP_T)cn, U);
				VECTOR_ID(free)(U);
				U = temp;
			}
			
			MATRIX_ID(free)(indg);
			MATRIX_ID(free)(Sit);
		} else {
			
			// U(1)=[];
			VECTOR_T* temp = VECTOR_ID(alloc)(U->size - 1);
			VECTOR_ID(view) U_subv = VECTOR_ID(subvector)(U, 1, U->size - 1);
			VECTOR_ID(memcpy)(temp, &U_subv.vector);
			VECTOR_ID(free)(U);
			U = temp;
		}
		
		MATRIX_ID(free)(S);
		
		// ind=find(Ci==U(1));
		VECTOR_T* _Ci_eq_U_0 = compare_elements(Ci, fp_equal, VECTOR_ID(get)(U, 0));
		VECTOR_ID(free)(ind);
		ind = find(_Ci_eq_U_0);
		VECTOR_ID(free)(_Ci_eq_U_0);
		
		if (ind != NULL) {
			
			// bg=B(ind,ind);
			MATRIX_T* bg = ordinal_index(B, ind, ind);
			
			// Bg=bg-diag(sum(bg));
			VECTOR_T* sum_bg = sum(bg);
			MATRIX_T* diag_sum_bg = diag(sum_bg);
			VECTOR_ID(free)(sum_bg);
			MATRIX_ID(free)(Bg);
			Bg = bg;
			MATRIX_ID(sub)(Bg, diag_sum_bg);
			MATRIX_ID(free)(diag_sum_bg);
			
			// Ng=length(ind);
			Ng = length(ind);
		}
	}
	
	VECTOR_ID(free)(U);
	VECTOR_ID(free)(ind);
	MATRIX_ID(free)(Bg);
	gsl_eigen_symmv_free(eig);
	
	return Ci;
}
