#include <cmath>

#include "bct.h"

/*
 * Counts occurrences of three-node structural motifs in a weighted graph.
 * Returns intensity and (optionally) coherence and motif counts.
 */
MATRIX_T* BCT_NAMESPACE::motif3struct_wei(const MATRIX_T* W, MATRIX_T** Q, MATRIX_T** F) {
	if (safe_mode) check_status(W, SQUARE | WEIGHTED, "motif3struct_wei");
	
	// load motif34lib M3 M3n ID3 N3
	VECTOR_T* ID3;
	VECTOR_T* N3;
	MATRIX_T* M3 = motif3generate(&ID3, &N3);
	
	// n=length(W);
	int n = length(W);
	
	// I=zeros(13,n);
	MATRIX_T* I = zeros(13, n);
	
	// Q=zeros(13,n);
	if (Q != NULL) {
		*Q = zeros(13, n);
	}
	
	// F=zeros(13,n);
	if (F != NULL) {
		*F = zeros(13, n);
	}
	
	// A=1*(W~=0);
	MATRIX_T* A = compare_elements(W, fp_not_equal, 0.0);
	
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
						
						// w=[W(v1,u) W(v2,u) W(u,v1) W(v2,v1) W(u,v2) W(v1,v2)];
						int WA_rows[] = { v1, v2, u, v2, u, v1 };
						int WA_cols[] = { u, u, v1, v1, v2, v2 };
						VECTOR_T* w = VECTOR_ID(alloc)(6);
						for (int i = 0; i < 6; i++) {
							VECTOR_ID(set)(w, i, MATRIX_ID(get)(W, WA_rows[i], WA_cols[i]));
						}
						
						// s=uint32(sum(10.^(5:-1:0).*[A(v1,u) A(v2,u) A(u,v1) A(v2,v1) A(u,v2) A(v1,v2)]));
						VECTOR_T* s = VECTOR_ID(alloc)(6);
						for (int i = 0; i < 6; i++) {
							VECTOR_ID(set)(s, i, MATRIX_ID(get)(A, WA_rows[i], WA_cols[i]));
						}
						
						// ind=(s==M3n);
						int ind = 0;
						for ( ; ind < (int)M3->size1; ind++) {
							VECTOR_ID(view) M3_row_ind = MATRIX_ID(row)(M3, ind);
							if (compare_vectors(s, &M3_row_ind.vector) == 0) {
								break;
							}
						}
						VECTOR_ID(free)(s);
						if (ind < (int)M3->size1) {
							
							// M=w.*M3(ind,:);
							VECTOR_T* M = VECTOR_ID(alloc)(M3->size2);
							MATRIX_ID(get_row)(M, M3, ind);
							VECTOR_ID(mul)(M, w);
							
							// id=ID3(ind);
							int id = (int)VECTOR_ID(get)(ID3, ind) - 1;
							
							// l=N3(ind);
							int l = (int)VECTOR_ID(get)(N3, ind);
							
							// x=sum(M,2)/l;
							FP_T x = sum(M) / (FP_T)l;
							
							// M(M==0)=1;
							VECTOR_T* M_eq_0 = compare_elements(M, fp_equal, 0.0);
							logical_index_assign(M, M_eq_0, 1.0);
							VECTOR_ID(free)(M_eq_0);
							
							// i=prod(M,2)^(1/l);
							FP_T i = std::pow(prod(M), (FP_T)1.0 / l);
							VECTOR_ID(free)(M);
							
							// q=i/x;
							FP_T q = i / x;
							
							// I(id,[u v1 v2])=I(id,[u v1 v2])+[i i i];
							// Q(id,[u v1 v2])=Q(id,[u v1 v2])+[q q q];
							// F(id,[u v1 v2])=F(id,[u v1 v2])+[1 1 1];
							int IQF_cols[] = { u, v1, v2 };
							for (int j = 0; j < 3; j++) {
								MATRIX_ID(set)(I, id, IQF_cols[j], MATRIX_ID(get)(I, id, IQF_cols[j]) + i);
								if (Q != NULL) {
									MATRIX_ID(set)(*Q, id, IQF_cols[j], MATRIX_ID(get)(*Q, id, IQF_cols[j]) + q);
								}
								if (F != NULL) {
									MATRIX_ID(set)(*F, id, IQF_cols[j], MATRIX_ID(get)(*F, id, IQF_cols[j]) + 1.0);
								}
							}
						}
						
						VECTOR_ID(free)(w);
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
	VECTOR_ID(free)(N3);
	MATRIX_ID(free)(M3);
	MATRIX_ID(free)(A);
	MATRIX_ID(free)(As);
	return I;
}

/*
 * Returns per-motif metrics instead of per-motif, per-node metrics.
 */
VECTOR_T* BCT_NAMESPACE::motif3struct_wei_v(const MATRIX_T* W, VECTOR_T** Q, VECTOR_T** F) {
	MATRIX_T* _Q;
	MATRIX_T* _F;
	MATRIX_T* _I = motif3struct_wei(W, &_Q, &_F);
	if (Q != NULL) {
		*Q = sum(_Q, 2);
	}
	MATRIX_ID(free)(_Q);
	if (F != NULL) {
		*F = sum(_F, 2);
	}
	MATRIX_ID(free)(_F);
	VECTOR_T* I = sum(_I, 2);
	MATRIX_ID(free)(_I);
	return I;
}
