#include "bct.h"

/*
 * Computes node betweenness for a binary graph.
 */
VECTOR_T* BCT_NAMESPACE::betweenness_bin(const MATRIX_T* G) {
	VECTOR_T* BC;
	MATRIX_T* EBC = edge_betweenness_bin(G, &BC);
	MATRIX_ID(free)(EBC);
	return BC;
}

/*
 * Computes node and edge betweenness for a binary graph.
 */
MATRIX_T* BCT_NAMESPACE::edge_betweenness_bin(const MATRIX_T* G, VECTOR_T** BC) {
	if (safe_mode) check_status(G, SQUARE | BINARY, "edge_betweenness_bin");
	
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
		
		// D=false(1,n); D(u) = 1;
		VECTOR_T* D = VECTOR_ID(calloc)(n);
		VECTOR_ID(set)(D, u, 1.0);
		
		// NP=zeros(1,n); NP(u)=1;
		VECTOR_T* NP = zeros_vector(n);
		VECTOR_ID(set)(NP, u, 1.0);
		
		// P=false(n);
		MATRIX_T* P = MATRIX_ID(calloc)(n, n);
		
		// Q=zeros(1,n); q=n;
		VECTOR_T* Q = zeros_vector(n);
		int q = n - 1;
		
		// Gu=G;
		MATRIX_T* Gu = copy(G);
		
		// V=u;
		VECTOR_T* V = VECTOR_ID(alloc)(1);
		VECTOR_ID(set)(V, 0, (FP_T)u);
		
		// while V
		while (V != NULL) {
			
			// Gu(:,V)=0;
			VECTOR_T* Gu_rows = sequence(0, Gu->size1 - 1);
			ordinal_index_assign(Gu, Gu_rows, V, 0.0);
			VECTOR_ID(free)(Gu_rows);
			
			// for v=V
			for (int i_V = 0; i_V < (int)V->size; i_V++) {
				int v = (int)VECTOR_ID(get)(V, i_V);
				
				// Q(q)=v; q=q-1;
				VECTOR_ID(set)(Q, q--, (FP_T)v);
				
				// W=find(Gu(v,:));
				VECTOR_ID(view) Gu_row_v = MATRIX_ID(row)(Gu, v);
				VECTOR_T* W = find(&Gu_row_v.vector);
				if (W != NULL) {
					
					// for w=W
					for (int i_W = 0; i_W < (int)W->size; i_W++) {
						int w = (int)VECTOR_ID(get)(W, i_W);
						
						// if D(w)
						if (fp_nonzero(VECTOR_ID(get)(D, w))) {
							
							// NP(w)=NP(w)+NP(v);
							VECTOR_ID(set)(NP, w, VECTOR_ID(get)(NP, w) + VECTOR_ID(get)(NP, v));
							
							// P(w,v)=1;
							MATRIX_ID(set)(P, w, v, 1.0);
							
							// else
						} else {
							
							// D(w)=1;
							VECTOR_ID(set)(D, w, 1.0);
							
							// NP(w)=NP(v);
							VECTOR_ID(set)(NP, w, VECTOR_ID(get)(NP, v));
							
							// P(w,v)=1;
							MATRIX_ID(set)(P, w, v, 1.0);
						}
					}
					VECTOR_ID(free)(W);
				}
			}
			
			// V=find(any(Gu(V,:),1));
			VECTOR_T* Gu_cols = sequence(0, G->size2 - 1);
			MATRIX_T* Gu_idx = ordinal_index(Gu, V, Gu_cols);
			VECTOR_ID(free)(Gu_cols);
			VECTOR_T* any_Gu_idx = any(Gu_idx);
			MATRIX_ID(free)(Gu_idx);
			VECTOR_ID(free)(V);
			V = find(any_Gu_idx);
			VECTOR_ID(free)(any_Gu_idx);
		}
		
		MATRIX_ID(free)(Gu);
		
		// if ~all(D)
		if (all(D) == 0) {
			
			// Q(1:q)=find(~D);
			VECTOR_T* not_D = logical_not(D);
			VECTOR_T* find_not_D = find(not_D);
			VECTOR_ID(free)(not_D);
			VECTOR_ID(view) Q_subv = VECTOR_ID(subvector)(Q, 0, q + 1);
			VECTOR_ID(memcpy)(&Q_subv.vector, find_not_D);
			VECTOR_ID(free)(find_not_D);
		}
		
		VECTOR_ID(free)(D);
		
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
