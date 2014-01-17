#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Computes node betweenness for a weighted graph.
 */
VECTOR_T* BCT_NAMESPACE::betweenness_wei(const MATRIX_T* G) {
	VECTOR_T* BC;
	MATRIX_T* EBC = edge_betweenness_wei(G, &BC);
	MATRIX_ID(free)(EBC);
	return BC;
}

/*
 * Computes node and edge betweenness for a weighted graph.
 */
MATRIX_T* BCT_NAMESPACE::edge_betweenness_wei(const MATRIX_T* G, VECTOR_T** BC) {
	if (safe_mode) check_status(G, SQUARE | WEIGHTED, "edge_betweenness_wei");
	
	// n=length(G);
	int n = length(G);
	
	// BC=zeros(n,1);
	if (BC != NULL) {
		*BC = zeros_vector(n);
	}
	
	// EBC=zeros(n);
	MATRIX_T* EBC = zeros(n);
	
	// for u=1:n
	for (int u = 0; u < n; u++) {
		
		// D=inf(1,n); D(u) = 0;
		VECTOR_T* D = VECTOR_ID(alloc)(n);
		VECTOR_ID(set_all)(D, GSL_POSINF);
		VECTOR_ID(set)(D, u, 0.0);
		
		// NP=zeros(1,n); NP(u)=1;
		VECTOR_T* NP = zeros_vector(n);
		VECTOR_ID(set)(NP, u, 1.0);
		
		// S=true(1,n);
		VECTOR_T* S = VECTOR_ID(alloc)(n);
		VECTOR_ID(set_all)(S, 1.0);
		
		// P=false(n);
		MATRIX_T* P = MATRIX_ID(calloc)(n, n);
		
		// Q=zeros(1,n); q=n;
		VECTOR_T* Q = zeros_vector(n);
		int q = n - 1;
		
		// G1=G;
		MATRIX_T* G1 = copy(G);
		
		// V=u;
		VECTOR_T* V = VECTOR_ID(alloc)(1);
		VECTOR_ID(set)(V, 0, (FP_T)u);
		
		// while 1
		while (true) {
			
			// S(V)=0;
			ordinal_index_assign(S, V, 0.0);
			
			// G1(:,V)=0;
			VECTOR_T* G1_rows = sequence(0, G1->size1 - 1);
			ordinal_index_assign(G1, G1_rows, V, 0.0);
			VECTOR_ID(free)(G1_rows);
			
			// for v=V
			for (int i_V = 0; i_V < (int)V->size; i_V++) {
				int v = (int)VECTOR_ID(get)(V, i_V);
				
				// Q(q)=v; q=q-1;
				VECTOR_ID(set)(Q, q--, (FP_T)v);
				
				// W=find(G1(v,:));
				VECTOR_ID(view) G1_row_v = MATRIX_ID(row)(G1, v);
				VECTOR_T* W = find(&G1_row_v.vector);
				if (W != NULL) {
					
					// for w=W
					for (int i_W = 0; i_W < (int)W->size; i_W++) {
						int w = (int)VECTOR_ID(get)(W, i_W);
						
						// Duw=D(v)+G1(v,w);
						FP_T Duw = VECTOR_ID(get)(D, v) + MATRIX_ID(get)(G1, v, w);
						
						// if Duw<D(w)
						if (Duw < VECTOR_ID(get)(D, w)) {
							
							// D(w)=Duw;
							VECTOR_ID(set)(D, w, Duw);
							
							// NP(w)=NP(v);
							VECTOR_ID(set)(NP, w, VECTOR_ID(get)(NP, v));
							
							// P(w,:)=0;
							VECTOR_ID(view) P_row_w = MATRIX_ID(row)(P, w);
							VECTOR_ID(set_zero)(&P_row_w.vector);
							
							// P(w,v)=1;
							MATRIX_ID(set)(P, w, v, 1.0);
							
							// elseif Duw==D(w)
						} else if (fp_equal(Duw, VECTOR_ID(get)(D, w))) {
							
							// NP(w)=NP(w)+NP(v);
							VECTOR_ID(set)(NP, w, VECTOR_ID(get)(NP, w) + VECTOR_ID(get)(NP, v));
							
							// P(w,v)=1;
							MATRIX_ID(set)(P, w, v, 1.0);
						}
					}
					VECTOR_ID(free)(W);
				}
			}
			
			// if isempty(minD), break
			if (nnz(S) == 0) {
				break;
			} else {
				
				// minD=min(D(S))
				VECTOR_T* D_S = logical_index(D, S);
				FP_T minD = VECTOR_ID(min)(D_S);
				VECTOR_ID(free)(D_S);
				
				// elseif isinf(minD),
				if (gsl_isinf(minD) == 1) {
					
					// Q(1:q)=find(isinf(D)); break
					VECTOR_T* isinf_D = compare_elements(D, fp_equal, GSL_POSINF);
					VECTOR_T* find_isinf_D = find(isinf_D);
					VECTOR_ID(free)(isinf_D);
					VECTOR_ID(view) Q_subv = VECTOR_ID(subvector)(Q, 0, q + 1);
					VECTOR_ID(memcpy)(&Q_subv.vector, find_isinf_D);
					VECTOR_ID(free)(find_isinf_D);
					break;
				}
				
				// V=find(D==minD);
				VECTOR_ID(free)(V);
				VECTOR_T* D_eq_minD = compare_elements(D, fp_equal, minD);
				V = find(D_eq_minD);
				VECTOR_ID(free)(D_eq_minD);
			}
		}
		
		VECTOR_ID(free)(D);
		VECTOR_ID(free)(S);
		MATRIX_ID(free)(G1);
		VECTOR_ID(free)(V);
		
		// DP=zeros(n,1);
		VECTOR_T* DP = zeros_vector(n);
		
		// for w=Q(1:n-1);
		for (int i_Q = 0; i_Q < n - 1; i_Q++) {
			int w = (int)VECTOR_ID(get)(Q, i_Q);
			
			// BC(w)=BC(w)+DP(w)
			if (BC != NULL) {
				VECTOR_ID(set)(*BC, w, VECTOR_ID(get)(*BC, w) + VECTOR_ID(get)(DP, w));
			}
			
			// for v=find(P(w,:))
			VECTOR_ID(view) P_row_w = MATRIX_ID(row)(P, w);
			VECTOR_T* find_P_row_w = find(&P_row_w.vector);
			if (find_P_row_w != NULL) {
				for (int i_find_P_row_w = 0; i_find_P_row_w < (int)find_P_row_w->size; i_find_P_row_w++) {
					int v = (int)VECTOR_ID(get)(find_P_row_w, i_find_P_row_w);
					
					// DPvw=(1+DP(w)).*NP(v)./NP(w);
					FP_T DP_w = VECTOR_ID(get)(DP, w);
					FP_T NP_v = VECTOR_ID(get)(NP, v);
					FP_T NP_w = VECTOR_ID(get)(NP, w);
					FP_T DPvw = (1 + DP_w) * NP_v / NP_w;
					
					// DP(v)=DP(v)+DPvw;
					VECTOR_ID(set)(DP, v, VECTOR_ID(get)(DP, v) + DPvw);
					
					// EBC(v,w)=EBC(v,w)+DPvw;
					MATRIX_ID(set)(EBC, v, w, MATRIX_ID(get)(EBC, v, w) + DPvw);
				}
				VECTOR_ID(free)(find_P_row_w);
			}
		}
		
		VECTOR_ID(free)(NP);
		MATRIX_ID(free)(P);
		VECTOR_ID(free)(Q);
		VECTOR_ID(free)(DP);
	}
	
	return EBC;
}
