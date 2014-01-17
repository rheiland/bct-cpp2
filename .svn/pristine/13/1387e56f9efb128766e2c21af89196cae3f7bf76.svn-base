#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Computes the distance matrix for a weighted graph.
 */
MATRIX_T* BCT_NAMESPACE::distance_wei(const MATRIX_T* G) {
	if (safe_mode) check_status(G, SQUARE | WEIGHTED, "distance_wei");
	
	// n=length(G);
	int n = length(G);
	
	// D=zeros(n); D(~eye(n))=inf;
	MATRIX_T* D = MATRIX_ID(alloc)(n, n);
	MATRIX_ID(set_all)(D, GSL_POSINF);
	VECTOR_ID(view) diag_D = MATRIX_ID(diagonal)(D);
	VECTOR_ID(set_all)(&diag_D.vector, 0.0);
	
	// for u=1:n
	for (int u = 0; u < n; u++) {
		
		// S=true(1,n);
		VECTOR_T* S = VECTOR_ID(alloc)(n);
		VECTOR_ID(set_all)(S, 1.0);
		
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
				
				// W=find(G1(v,:));
				VECTOR_ID(view) G1_row_v = MATRIX_ID(row)(G1, v);
				VECTOR_T* W = find(&G1_row_v.vector);
				if (W != NULL) {
					
					// D(u,W)=min([D(u,W);D(u,v)+G1(v,W)]);
					VECTOR_T* DG1_row = VECTOR_ID(alloc)(1);
					VECTOR_ID(set)(DG1_row, 0, (FP_T)u);
					MATRIX_T* D_u_W = ordinal_index(D, DG1_row, W);
					FP_T D_u_v = MATRIX_ID(get)(D, u, v);
					VECTOR_ID(set)(DG1_row, 0, (FP_T)v);
					MATRIX_T* G1_v_W = ordinal_index(G1, DG1_row, W);
					MATRIX_ID(add_constant)(G1_v_W, D_u_v);
					MATRIX_T* cat = concatenate_columns(D_u_W, G1_v_W);
					MATRIX_ID(free)(D_u_W);
					MATRIX_ID(free)(G1_v_W);
					VECTOR_T* min_cat_v = min(cat);
					MATRIX_ID(free)(cat);
					MATRIX_T* min_cat = to_row_matrix(min_cat_v);
					VECTOR_ID(free)(min_cat_v);
					VECTOR_ID(set)(DG1_row, 0, (FP_T)u);
					ordinal_index_assign(D, DG1_row, W, min_cat);
					VECTOR_ID(free)(W);
					VECTOR_ID(free)(DG1_row);
					MATRIX_ID(free)(min_cat);
				}
			}
			
			// minD=min(D(u,S));
			// if isempty(minD)||isinf(minD), break, end;
			VECTOR_T* D_row = VECTOR_ID(alloc)(1);
			VECTOR_ID(set)(D_row, 0, (FP_T)u);
			MATRIX_T* D_u_S = ord_log_index(D, D_row, S);
			VECTOR_ID(free)(D_row);
			if (D_u_S == NULL) {
				break;
			} else {
				FP_T minD = MATRIX_ID(min)(D_u_S);
				MATRIX_ID(free)(D_u_S);
				if (gsl_isinf(minD) == 1) {
					break;
				}
				
				// V=find(D(u,:)==minD);
				VECTOR_ID(free)(V);
				VECTOR_ID(view) D_row_u = MATRIX_ID(row)(D, u);
				VECTOR_T* D_row_u_eq_minD = compare_elements(&D_row_u.vector, fp_equal, minD);
				V = find(D_row_u_eq_minD);
				VECTOR_ID(free)(D_row_u_eq_minD);
			}
		}
		
		VECTOR_ID(free)(S);
		MATRIX_ID(free)(G1);
		VECTOR_ID(free)(V);
	}
	
	return D;
}
