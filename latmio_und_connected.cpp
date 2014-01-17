#include <cmath>

#include "bct.h"

/*
 * Returns a latticized graph with equivalent degree sequence to the original
 * weighted undirected graph, and with preserved connectedness.  On average,
 * each edge is rewired ITER times.  Strength distributions are not preserved
 * for weighted graphs.
 */
MATRIX_T* BCT_NAMESPACE::latmio_und_connected(const MATRIX_T* R, int ITER) {
	if (safe_mode) check_status(R, SQUARE | UNDIRECTED, "latmio_und_connected");
	
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
	
	// [i j]=find(tril(R));
	MATRIX_T* tril_R = tril(R);
	MATRIX_T* find_tril_R = find_ij(tril_R);
	MATRIX_ID(free)(tril_R);
	VECTOR_ID(view) i = MATRIX_ID(column)(find_tril_R, 0);
	VECTOR_ID(view) j = MATRIX_ID(column)(find_tril_R, 1);
	
	// K=length(i);
	int K = length(&i.vector);
	
	// ITER=K*ITER;
	ITER = K * ITER;
	
	MATRIX_T* _R = copy(R);
	
	// for iter=1:ITER
	for (int iter = 1; iter <= ITER; iter++) {
		
		// while 1
		while (true) {
			
			// rewire = 1
			bool rewire = true;
			
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
			
			// if rand>0.5
			if (gsl_rng_uniform(rng) > 0.5) {
				
				// i(e2)=d; j(e2)=c;
				VECTOR_ID(set)(&i.vector, e2, (FP_T)d);
				VECTOR_ID(set)(&j.vector, e2, (FP_T)c);
				
				// c=i(e2); d=j(e2);
				c = (int)VECTOR_ID(get)(&i.vector, e2);
				d = (int)VECTOR_ID(get)(&j.vector, e2);
			}
			
			// if ~(R(a,d) || R(c,b))
			if (fp_zero(MATRIX_ID(get)(_R, a, d)) && fp_zero(MATRIX_ID(get)(_R, c, b))) {
				
				// if (D(a,b)+D(c,d))>=(D(a,d)+D(c,b))
				if (fp_greater_or_equal(MATRIX_ID(get)(D, a, b) + MATRIX_ID(get)(D, c, d),
										MATRIX_ID(get)(D, a, d) + MATRIX_ID(get)(D, c, b))) {
					
					// if ~(R(a,c) || R(b,d))
					if (fp_zero(MATRIX_ID(get)(_R, a, c)) && fp_zero(MATRIX_ID(get)(_R, b, d))) {
						
						// P=R([a d],:);
						VECTOR_T* _R_rows = VECTOR_ID(alloc)(2);
						VECTOR_ID(set)(_R_rows, 0, (FP_T)a);
						VECTOR_ID(set)(_R_rows, 1, (FP_T)d);
						VECTOR_T* _R_cols = sequence(0, _R->size2 - 1);
						MATRIX_T* P = ordinal_index(_R, _R_rows, _R_cols);
						VECTOR_ID(free)(_R_rows);
						VECTOR_ID(free)(_R_cols);
						
						// P(1,b)=0; P(2,c)=0;
						MATRIX_ID(set)(P, 0, b, 0.0);
						MATRIX_ID(set)(P, 1, c, 0.0);
						
						// PN=P;
						MATRIX_T* PN = copy(P);
						
						// PN(:,d)=1; PN(:,a)=1;
						VECTOR_ID(view) PN_col_d = MATRIX_ID(column)(PN, d);
						VECTOR_ID(set_all)(&PN_col_d.vector, 1.0);
						VECTOR_ID(view) PN_col_a = MATRIX_ID(column)(PN, a);
						VECTOR_ID(set_all)(&PN_col_a.vector, 1.0);
						
						// while 1
						while (true) {
							
							// P(1,:)=any(R(P(1,:)~=0,:),1);
							VECTOR_ID(view) P_row_0 = MATRIX_ID(row)(P, 0);
							VECTOR_T* P_row_0_neq_0 = compare_elements(&P_row_0.vector, fp_not_equal, 0.0);
							VECTOR_T* _R_cols = sequence(0, _R->size2 - 1);
							MATRIX_T* _R_idx = log_ord_index(_R, P_row_0_neq_0, _R_cols);
							VECTOR_ID(free)(P_row_0_neq_0);
							if (_R_idx != NULL) {
								VECTOR_T* any__R_idx = any(_R_idx, 1);
								MATRIX_ID(free)(_R_idx);
								MATRIX_ID(set_row)(P, 0, any__R_idx);
								VECTOR_ID(free)(any__R_idx);
							} else {
								VECTOR_ID(set_zero)(&P_row_0.vector);
							}
							
							// P(2,:)=any(R(P(2,:)~=0,:),1);
							VECTOR_ID(view) P_row_1 = MATRIX_ID(row)(P, 0);
							VECTOR_T* P_row_1_neq_0 = compare_elements(&P_row_1.vector, fp_not_equal, 0.0);
							_R_idx = log_ord_index(_R, P_row_1_neq_0, _R_cols);
							VECTOR_ID(free)(P_row_1_neq_0);
							VECTOR_ID(free)(_R_cols);
							if (_R_idx != NULL) {
								VECTOR_T* any__R_idx = any(_R_idx, 1);
								MATRIX_ID(free)(_R_idx);
								MATRIX_ID(set_row)(P, 1, any__R_idx);
								VECTOR_ID(free)(any__R_idx);
							} else {
								VECTOR_ID(set_zero)(&P_row_1.vector);
							}
							
							// P=P.*(~PN);
							MATRIX_T* not_PN = logical_not(PN);
							MATRIX_ID(mul_elements)(P, not_PN);
							MATRIX_ID(free)(not_PN);
							
							// if ~all(any(P,2))
							VECTOR_T* any_P = any(P, 2);
							bool all_any_P = all(any_P);
							VECTOR_ID(free)(any_P);
							if (!all_any_P) {
								
								// rewire=0;
								rewire = false;
								
								// break
								break;
							}
							
							// elseif any(any(P(:,[b c])))
							VECTOR_ID(view) P_col_b = MATRIX_ID(column)(P, b);
							VECTOR_ID(view) P_col_c = MATRIX_ID(column)(P, c);
							MATRIX_T* P_idx = concatenate_rows(&P_col_b.vector, &P_col_c.vector);
							VECTOR_T* any_P_idx = any(P_idx);
							MATRIX_ID(free)(P_idx);
							bool any_any_P_idx = any(any_P_idx);
							VECTOR_ID(free)(any_P_idx);
							if (any_any_P_idx) {
								
								// break
								break;
							}
							
							// PN=PN+P;
							MATRIX_ID(add)(PN, P);
						}
						
						MATRIX_ID(free)(P);
						MATRIX_ID(free)(PN);
					}
					
					// if rewire
					if (rewire) {
						
						// R(a,d)=R(a,b); R(a,b)=0;
						MATRIX_ID(set)(_R, a, d, MATRIX_ID(get)(_R, a, b));
						MATRIX_ID(set)(_R, a, b, 0.0);
						
						// R(d,a)=R(b,a); R(b,a)=0;
						MATRIX_ID(set)(_R, d, a, MATRIX_ID(get)(_R, b, a));
						MATRIX_ID(set)(_R, b, a, 0.0);
						
						// R(c,b)=R(c,d); R(c,d)=0;
						MATRIX_ID(set)(_R, c, b, MATRIX_ID(get)(_R, c, d));
						MATRIX_ID(set)(_R, c, d, 0.0);
						
						// R(b,c)=R(d,c); R(d,c)=0;
						MATRIX_ID(set)(_R, b, c, MATRIX_ID(get)(_R, d, c));
						MATRIX_ID(set)(_R, d, c, 0.0);
						
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
	}
	
	MATRIX_ID(free)(D);
	MATRIX_ID(free)(find_tril_R);
	return _R;
}
