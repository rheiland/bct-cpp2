#include "bct.h"

bool modularity_louvain_und(const MATRIX_T*, FP_T*, VECTOR_T**, int);

/*
 * Detects communities in an undirected graph via Louvain modularity.  While the
 * MATLAB version returns intermediate values for community numbering and
 * modularity, the C++ version returns only the final values.  This version also
 * makes use of an additional argument N that specifies the maximum number of
 * node permutations to attempt when maximizing modularity.
 */
FP_T BCT_NAMESPACE::modularity_louvain_und(const MATRIX_T* W, VECTOR_T** Ci, int N) {
	if (safe_mode) check_status(W, SQUARE | UNDIRECTED, "modularity_louvain_und");
	
	FP_T Q;
	while (true) {
		if (modularity_louvain_und(W, &Q, Ci, N)) {
			break;
		}
	}
	return Q;
}

bool modularity_louvain_und(const MATRIX_T* W, FP_T* Q, VECTOR_T** Ci, int N) {
	using namespace BCT_NAMESPACE;
	
	// n=length(W);
	int n = length(W);
	
	// s=sum(W(:));
	VECTOR_T* sum_W = sum(W);
	FP_T s = sum(sum_W);
	VECTOR_ID(free)(sum_W);
	
	// h=1;
	int h = 0;
	
	// Ci{h}=1:n;
	std::vector<VECTOR_T*> _Ci;
	_Ci.push_back(sequence(1, n));
	
	// Q{h}=-1;
	std::vector<FP_T> _Q;
	_Q.push_back(-1.0);
	
	// n0=n;
	int n0 = n;
	
	MATRIX_T* _W = copy(W);
	
	// while true
	while (true) {
		
		// K=sum(W);
		VECTOR_T* K = sum(_W);
		
		// Km=K;
		VECTOR_T* Km = copy(K);
		
		// Knm=W;
		MATRIX_T* Knm = copy(_W);
		
		// M=1:n;
		VECTOR_T* M = sequence(0, n - 1);
		
		// Nm=ones(1,n);
		VECTOR_T* Nm = ones_vector(n);
		
		// flag=true;
		bool flag = true;
		
		int count = 0;
		
		// while flag
		while (flag) {
			
			if (++count >= N) {
				return false;
			}

			// flag=false;
			flag = false;
			
			// for i=randperm(n)
			gsl_permutation* randperm_n = randperm(n);
			for (int i_randperm_n = 0; i_randperm_n < n; i_randperm_n++) {
				int i = gsl_permutation_get(randperm_n, i_randperm_n);
				
				// dQ=(Knm(i,:)-Knm(i,M(i))+W(i,i)) - K(i).*(Km-Km(M(i))+K(i))/s;
				VECTOR_T* dQ1 = VECTOR_ID(alloc)(Knm->size2);
				MATRIX_ID(get_row)(dQ1, Knm, i);
				int M_i = (int)VECTOR_ID(get)(M, i);
				VECTOR_ID(add_constant)(dQ1, -MATRIX_ID(get)(Knm, i, M_i));
				VECTOR_ID(add_constant)(dQ1, MATRIX_ID(get)(_W, i, i));
				VECTOR_T* dQ2 = copy(Km);
				VECTOR_ID(add_constant)(dQ2, -VECTOR_ID(get)(Km, M_i));
				VECTOR_ID(add_constant)(dQ2, VECTOR_ID(get)(K, i));
				VECTOR_ID(scale)(dQ2, VECTOR_ID(get)(K, i) / s);
				VECTOR_T* dQ = dQ1;
				VECTOR_ID(sub)(dQ, dQ2);
				VECTOR_ID(free)(dQ2);
				
				// dQ(M(i))=0;
				VECTOR_ID(set)(dQ, M_i, 0.0);
				
				// max_dQ=max(dQ);
				FP_T max_dQ = max(dQ);
				
				// if max_dQ>0;
				if (max_dQ > 0.0) {
					
					// j=find(dQ==max_dQ,1);
					VECTOR_T* dQ_eq_max_dQ = compare_elements(dQ, fp_equal, max_dQ);
					VECTOR_T* j_v = find(dQ_eq_max_dQ, 1);
					VECTOR_ID(free)(dQ_eq_max_dQ);
					int j = (int)VECTOR_ID(get)(j_v, 0);
					VECTOR_ID(free)(j_v);
					
					// Knm(:,j)=Knm(:,j)+W(:,i);
					VECTOR_ID(view) Knm_col_j = MATRIX_ID(column)(Knm, j);
					VECTOR_ID(view) _W_col_i = MATRIX_ID(column)(_W, i);
					VECTOR_ID(add)(&Knm_col_j.vector, &_W_col_i.vector);
					
					// Knm(:,M(i))=Knm(:,M(i))-W(:,i);
					VECTOR_ID(view) Knm_col_M_i = MATRIX_ID(column)(Knm, M_i);
					VECTOR_ID(sub)(&Knm_col_M_i.vector, &_W_col_i.vector);
					
					// Km(j)=Km(j)+K(i);
					VECTOR_ID(set)(Km, j, VECTOR_ID(get)(Km, j) + VECTOR_ID(get)(K, i));
					
					// Km(M(i))=Km(M(i))-K(i);
					VECTOR_ID(set)(Km, M_i, VECTOR_ID(get)(Km, M_i) - VECTOR_ID(get)(K, i));
					
					// Nm(j)=Nm(j)+1;
					VECTOR_ID(set)(Nm, j, VECTOR_ID(get)(Nm, j) + 1.0);
					
					// Nm(M(i))=Nm(M(i))-1;
					VECTOR_ID(set)(Nm, M_i, VECTOR_ID(get)(Nm, M_i) - 1.0);
					
					// M(i)=j;
					VECTOR_ID(set)(M, i, (FP_T)j);
					
					// flag=true;
					flag = true;
				}
				
				VECTOR_ID(free)(dQ);
			}
			
			gsl_permutation_free(randperm_n);
		}
		
		VECTOR_ID(free)(K);
		VECTOR_ID(free)(Km);
		MATRIX_ID(free)(Knm);
		VECTOR_ID(free)(Nm);
		
		// [x x M1]=unique(M);
		VECTOR_T* x1;
		VECTOR_T* M1;
		VECTOR_T* x2 = unique(M, "last", &x1, &M1);
		VECTOR_ID(free)(M);
		VECTOR_ID(free)(x1);
		VECTOR_ID(free)(x2);
		
		// h=h+1;
		h++;
		
		// Ci{h}=zeros(1,n0);
		_Ci.push_back(zeros_vector(n0));
		
		// for i=1:n
		for (int i = 0; i < n; i++) {
			
			// Ci{h}(Ci{h-1}==i)=M1(i);
			VECTOR_T* _Ci_h_sub_1_eq_i_add_1 = compare_elements(_Ci[h - 1], fp_equal, (FP_T)(i + 1));
			logical_index_assign(_Ci[h], _Ci_h_sub_1_eq_i_add_1, VECTOR_ID(get)(M1, i) + 1.0);
			VECTOR_ID(free)(_Ci_h_sub_1_eq_i_add_1);
		}
		
		// n=max(M1);
		n = (int)max(M1) + 1;
		
		// W1=zeros(n);
		MATRIX_T* _W1 = MATRIX_ID(alloc)(n, n);
		
		// for i=1:n
		for (int i = 0; i < n; i++) {
			
			// for j=1:n
			for (int j = 0; j < n; j++) {
				
				// w=sum(sum(W(M1==i,M1==j)));
				VECTOR_T* M1_eq_i = compare_elements(M1, fp_equal, (FP_T)i);
				VECTOR_T* M1_eq_j = compare_elements(M1, fp_equal, (FP_T)j);
				MATRIX_T* _W_idx = logical_index(_W, M1_eq_i, M1_eq_j);
				VECTOR_ID(free)(M1_eq_i);
				VECTOR_ID(free)(M1_eq_j);
				VECTOR_T* sum__W_idx = sum(_W_idx);
				MATRIX_ID(free)(_W_idx);
				FP_T w = sum(sum__W_idx);
				VECTOR_ID(free)(sum__W_idx);
				
				// W1(i,j)=w;
				MATRIX_ID(set)(_W1, i, j, w);
				
				// W1(j,i)=w;
				MATRIX_ID(set)(_W1, j, i, w);
			}
		}
		
		VECTOR_ID(free)(M1);
		
		// W=W1;
		MATRIX_ID(free)(_W);
		_W = _W1;
		
		// Q{h}=sum(diag(W))/s-sum(sum((W/s)^2));
		VECTOR_T* diag__W = diag(_W);
		FP_T sum_diag__W = sum(diag__W);
		VECTOR_ID(free)(diag__W);
		MATRIX_T* _W_div_s = copy(_W);
		MATRIX_ID(scale)(_W_div_s, 1.0 / s);
		MATRIX_T* _W_div_s_pow_2 = pow(_W_div_s, 2);
		MATRIX_ID(free)(_W_div_s);
		VECTOR_T* sum__W_div_s_pow_2 = sum(_W_div_s_pow_2);
		MATRIX_ID(free)(_W_div_s_pow_2);
		FP_T sum_sum__W_div_s_pow_2 = sum(sum__W_div_s_pow_2);
		VECTOR_ID(free)(sum__W_div_s_pow_2);
		_Q.push_back(sum_diag__W / s - sum_sum__W_div_s_pow_2);
		
		// if Q{h}-Q{h-1}<=eps
		if (fp_less_or_equal(_Q[h] - _Q[h - 1], epsilon)) {
			
			// break
			break;
		}
	}
				
	MATRIX_ID(free)(_W);
	
	// Ci([1 end])=[];
	for (int i = 0; i < (int)_Ci.size(); i++) {
		if (i != (int)_Ci.size() - 2) {
			VECTOR_ID(free)(_Ci[i]);
		}
	}
	if (Ci == NULL) {
		VECTOR_ID(free)(_Ci[_Ci.size() - 2]);
	} else {
		*Ci = _Ci[_Ci.size() - 2];
	} 
	
	// Q([1 end])=[];
	*Q = _Q[_Q.size() - 2];
	
	return true;
}
