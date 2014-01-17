#include "bct.h"

/*
 * Counts occurrences of three-node structural motifs in a binary graph.
 */
VECTOR_T* BCT_NAMESPACE::motif3struct_bin(const MATRIX_T* A, MATRIX_T** F) {
	if (safe_mode) check_status(A, SQUARE | BINARY, "motif3struct_bin");
	
	// load motif34lib M3n ID3
	VECTOR_T* ID3;
	MATRIX_T* M3 = motif3generate(&ID3);
	
	// n=length(A);
	int n = length(A);
	
	// F=zeros(13,n);
	if (F != NULL) {
		*F = zeros(13, n);
	}
	
	// f=zeros(13,1);
	VECTOR_T* f = zeros_vector(13);
	
	// As=A|A.';
	MATRIX_T* A_transpose = MATRIX_ID(alloc)(A->size2, A->size1);
	MATRIX_ID(transpose_memcpy)(A_transpose, A);
	MATRIX_T* As = logical_or(A, A_transpose);
	MATRIX_ID(free)(A_transpose);
	
	// for u=1:n-2
	for (int u = 0; u < n - 2; u++) {
		
		// V1=[false(1,u) As(u,u+1:n)];
		VECTOR_T* V1 = VECTOR_ID(alloc)(n);
		MATRIX_ID(get_row)(V1, As, u);
		for (int i = 0; i <= u; i++) {
			VECTOR_ID(set)(V1, i, 0.0);
		}
		
		// for v1=find(V1)
		VECTOR_T* find_V1 = find(V1);
		if (find_V1 != NULL) {
			for (int i_find_V1 = 0; i_find_V1 < (int)find_V1->size; i_find_V1++) {
				int v1 = (int)VECTOR_ID(get)(find_V1, i_find_V1);
				
				// V2=[false(1,u) As(v1,u+1:n)];
				VECTOR_T* V2 = VECTOR_ID(alloc)(n);
				MATRIX_ID(get_row)(V2, As, v1);
				for (int i = 0; i <= u; i++) {
					VECTOR_ID(set)(V2, i, 0.0);
				}
				
				// V2(V1)=0;
				logical_index_assign(V2, V1, 0.0);
				
				// V2=([false(1,v1) As(u,v1+1:n)])|V2;
				VECTOR_T* V2_1 = VECTOR_ID(alloc)(n);
				MATRIX_ID(get_row)(V2_1, As, u);
				for (int i = 0; i <= v1; i++) {
					VECTOR_ID(set)(V2_1, i, 0.0);
				}
				VECTOR_T* V2_2 = V2;
				V2 = logical_or(V2_1, V2_2);
				VECTOR_ID(free)(V2_1);
				VECTOR_ID(free)(V2_2);
				
				// for v2=find(V2)
				VECTOR_T* find_V2 = find(V2);
				if (find_V2 != NULL) {
					for (int i_find_V2 = 0; i_find_V2 < (int)find_V2->size; i_find_V2++) {
						int v2 = (int)VECTOR_ID(get)(find_V2, i_find_V2);
						
						// s=uint32(sum(10.^(5:-1:0).*[A(v1,u) A(v2,u) A(u,v1) A(v2,v1) A(u,v2) A(v1,v2)]));
						int A_rows[] = { v1, v2, u, v2, u, v1 };
						int A_cols[] = { u, u, v1, v1, v2, v2 };
						VECTOR_T* s = VECTOR_ID(alloc)(6);
						for (int i = 0; i < 6; i++) {
							VECTOR_ID(set)(s, i, MATRIX_ID(get)(A, A_rows[i], A_cols[i]));
						}
						
						// ind=ID3(s==M3n);
						int i_M3 = 0;
						for ( ; i_M3 < (int)M3->size1; i_M3++) {
							VECTOR_ID(view) M3_row_i_M3 = MATRIX_ID(row)(M3, i_M3);
							if (compare_vectors(s, &M3_row_i_M3.vector) == 0) {
								break;
							}
						}
						VECTOR_ID(free)(s);
						if (i_M3 < (int)M3->size1) {
							int ind = (int)VECTOR_ID(get)(ID3, i_M3) - 1;
							
							// if nargout==2; F(ind,[u v1 v2])=F(ind,[u v1 v2])+1; end
							if (F != NULL) {
								int F_cols[] = { u, v1, v2 };
								for (int i = 0; i < 3; i++) {
									MATRIX_ID(set)(*F, ind, F_cols[i], MATRIX_ID(get)(*F, ind, F_cols[i]) + 1.0);
								}
							}
							
							// f(ind)=f(ind)+1;
							VECTOR_ID(set)(f, ind, VECTOR_ID(get)(f, ind) + 1.0);
						}
					}
					
					VECTOR_ID(free)(find_V2);
				}
				
				VECTOR_ID(free)(V2);
			}
			
			VECTOR_ID(free)(find_V1);
		}
		
		VECTOR_ID(free)(V1);
	}
	
	VECTOR_ID(free)(ID3);
	MATRIX_ID(free)(M3);
	MATRIX_ID(free)(As);
	return f;
}
