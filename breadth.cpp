#include <gsl/gsl_math.h>

#include "bct.h"

/*
 * Performs a breadth-first search starting at the source node.  Because C++
 * indexing is zero-based, a value of 0 at branch(i) could mean either that node
 * 0 precedes node i or that node i is unreachable.  Check distance(i) for
 * GSL_POSINF to differentiate between these two cases.
 */
VECTOR_T* BCT_NAMESPACE::breadth(const MATRIX_T* CIJ, int source, VECTOR_T** branch) {
	if (safe_mode) check_status(CIJ, SQUARE, "breadth");
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// % colors: white, gray, black
	int white = 0;
	int gray = 1;
	int black = 2;
	
	// color = zeros(1,N);
	VECTOR_T* color = zeros_vector(N);
	
	// distance = inf*ones(1,N);
	VECTOR_T* distance = VECTOR_ID(alloc)(N);
	VECTOR_ID(set_all)(distance, GSL_POSINF);
	
	// branch = zeros(1,N);
	if (branch != NULL) {
		*branch = zeros_vector(N);
	}
	
	// color(source) = gray;
	VECTOR_ID(set)(color, source, (FP_T)gray);
	
	// distance(source) = 0;
	VECTOR_ID(set)(distance, source, 0.0);
	
	// branch(source) = -1;
	if (branch != NULL) {
		VECTOR_ID(set)(*branch, source, -1.0);
	}
	
	// Q = source;
	VECTOR_T* Q = VECTOR_ID(alloc)(1);
	VECTOR_ID(set)(Q, 0, (FP_T)source);
	
	// while ~isempty(Q)
	while (Q != NULL) {
		
		// u = Q(1);
		int u = (int)VECTOR_ID(get)(Q, 0);
		
		// ns = find(CIJ(u,:));
		VECTOR_ID(const_view) CIJ_row_u = MATRIX_ID(const_row)(CIJ, u);
		VECTOR_T* ns = find(&CIJ_row_u.vector);
		
		// for v=ns
		if (ns != NULL) {
			for (int i_ns = 0; i_ns < (int)ns->size; i_ns++) {
				int v = (int)VECTOR_ID(get)(ns, i_ns);
				
				// if (distance(v)==0)
				if ((int)VECTOR_ID(get)(distance, v) == 0) {
					
					// distance(v) = distance(u)+1;
					VECTOR_ID(set)(distance, v, VECTOR_ID(get)(distance, u) + 1.0);
				}
				
				// if (color(v)==white)
				if ((int)VECTOR_ID(get)(color, v) == white) {
					
					// color(v) = gray;
					VECTOR_ID(set)(color, v, (FP_T)gray);
					
					// distance(v) = distance(u)+1;
					VECTOR_ID(set)(distance, v, VECTOR_ID(get)(distance, u) + 1.0);
					
					// branch(v) = u;
					if (branch != NULL) {
						VECTOR_ID(set)(*branch, v, (FP_T)u);
					}
					
					// Q = [Q v];
					VECTOR_T* temp = concatenate(Q, (FP_T)v);
					VECTOR_ID(free)(Q);
					Q = temp;
				}
			}
			
			VECTOR_ID(free)(ns);
		}
		
		// Q = Q(2:length(Q));
		VECTOR_T* Q_cols = sequence(1, length(Q) - 1);
		if (Q_cols != NULL) {
			VECTOR_T* temp = ordinal_index(Q, Q_cols);
			VECTOR_ID(free)(Q);
			VECTOR_ID(free)(Q_cols);
			Q = temp;
		} else {
			VECTOR_ID(free)(Q);
			Q = NULL;
		}
		
		// color(u) = black;
		VECTOR_ID(set)(color, u, (FP_T)black);
	}
	
	VECTOR_ID(free)(color);
	return distance;
}
