#include <cmath>

#include "bct.h"

/*
 * Returns a latticized graph with equivalent degree sequence to the original
 * weighted directed graph.  On average, each edge is rewired ITER times.  Out-
 * strength is preserved for weighted graphs, while in-strength is not.
 */
MATRIX_T* BCT_NAMESPACE::latmio_dir(const MATRIX_T* R, int ITER) {
	if (safe_mode) check_status(R, SQUARE | DIRECTED, "latmio_dir");
	
	gsl_rng* rng = get_rng();
	
	// n=length(R);
	int n = length(R);
	
	// D=zeros(n);
	MATRIX_T* D = zeros(n);
	
	// u=[0 min([mod(1:n-1,n);mod(n-1:-1:1,n)])];
	VECTOR_T* seq1 = sequence(1, n - 1);
	VECTOR_T* seq2 = sequence(n - 1, -1, 1);
	MATRIX_T* seq1_seq2 = concatenate_columns(seq1, seq2);
	VECTOR_ID(free)(seq1);
	VECTOR_ID(free)(seq2);
	VECTOR_T* min_seq1_seq2 = min(seq1_seq2);
	MATRIX_ID(free)(seq1_seq2);
	VECTOR_T* u = concatenate(0.0, min_seq1_seq2);
	VECTOR_ID(free)(min_seq1_seq2);
	
	// for v=1:ceil(n/2)
	for (int v = 1; v <= (int)std::ceil((FP_T)n / 2.0); v++) {
		
		// D(n-v+1,:)=u([v+1:n 1:v]);
		VECTOR_T* u_indices1 = sequence(v, n - 1);
		VECTOR_T* u_indices2 = sequence(0, v - 1);
		VECTOR_T* u_indices = concatenate(u_indices1, u_indices2);
		VECTOR_ID(free)(u_indices1);
		VECTOR_ID(free)(u_indices2);
		VECTOR_T* u_idx = ordinal_index(u, u_indices);
		VECTOR_ID(free)(u_indices);
		MATRIX_ID(set_row)(D, n - v, u_idx);
		VECTOR_ID(free)(u_idx);
		
		// D(v,:)=D(n-v+1,n:-1:1);
		VECTOR_T* D_rows = VECTOR_ID(alloc)(1);
		VECTOR_ID(set)(D_rows, 0, (FP_T)(n - v));
		VECTOR_T* D_cols = sequence(n - 1, -1, 0);
		MATRIX_T* D_idx = ordinal_index(D, D_rows, D_cols);
		VECTOR_ID(free)(D_rows);
		VECTOR_ID(free)(D_cols);
		VECTOR_T* D_idx_v = to_vector(D_idx);
		MATRIX_ID(free)(D_idx);
		MATRIX_ID(set_row)(D, v - 1, D_idx_v);
		VECTOR_ID(free)(D_idx_v);
	}
	
	// [i j]=find(R);
	MATRIX_T* find_R = find_ij(R);
	VECTOR_ID(view) i = MATRIX_ID(column)(find_R, 0);
	VECTOR_ID(view) j = MATRIX_ID(column)(find_R, 1);
	
	// K=length(i);
	int K = length(&i.vector);
	
	// ITER=K*ITER;
	ITER = K * ITER;
	
	MATRIX_T* _R = copy(R);
	
	// for iter=1:ITER
	for (int iter = 1; iter <= ITER; iter++) {
		
		// while 1
		while (true) {
			
			int e1, e2;
			int a, b, c, d;
			
			// while 1
			while (true) {
				
				// e1=ceil(K*rand);
				e1 = gsl_rng_uniform_int(rng, K);
				
				// e2=ceil(K*rand);
				e2 = gsl_rng_uniform_int(rng, K);
				
				// while (e2==e1),
				while (e2 == e1) {
					
					// e2=ceil(K*rand);
					e2 = gsl_rng_uniform_int(rng, K);
				}
				
				// a=i(e1); b=j(e1);
				a = (int)VECTOR_ID(get)(&i.vector, e1);
				b = (int)VECTOR_ID(get)(&j.vector, e1);
				
				// c=i(e2); d=j(e2);
				c = (int)VECTOR_ID(get)(&i.vector, e2);
				d = (int)VECTOR_ID(get)(&j.vector, e2);
				
				// if all(a~=[c d]) && all(b~=[c d]);
				if (a != c && a != d && b != c && b != d) {
					
					// break
					break;
				}
			}
			
			// if ~(R(a,d) || R(c,b))
			if (fp_zero(MATRIX_ID(get)(_R, a, d)) && fp_zero(MATRIX_ID(get)(_R, c, b))) {
				
				// if (D(a,b)+D(c,d))>=(D(a,d)+D(c,b))
				if (fp_greater_or_equal(MATRIX_ID(get)(D, a, b) + MATRIX_ID(get)(D, c, d),
										MATRIX_ID(get)(D, a, d) + MATRIX_ID(get)(D, c, b))) {
					
					// R(a,d)=R(a,b); R(a,b)=0;
					MATRIX_ID(set)(_R, a, d, MATRIX_ID(get)(_R, a, b));
					MATRIX_ID(set)(_R, a, b, 0.0);
					
					// R(c,b)=R(c,d); R(c,d)=0;
					MATRIX_ID(set)(_R, c, b, MATRIX_ID(get)(_R, c, d));
					MATRIX_ID(set)(_R, c, d, 0.0);
					
					// j(e1) = d;
					VECTOR_ID(set)(&j.vector, e1, (FP_T)d);
					
					// j(e2) = b;
					VECTOR_ID(set)(&j.vector, e2, (FP_T)b);
					
					// break;
					break;
				}
			}
		}
	}
	
	MATRIX_ID(free)(D);
	MATRIX_ID(free)(find_R);
	return _R;
}
