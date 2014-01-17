#include <cmath>

#include "bct.h"

/*
 * Counts occurrences of four-node structural motifs in a weighted graph.
 * Returns intensity and (optionally) coherence and motif counts.
 */
MATRIX_T* BCT_NAMESPACE::motif4struct_wei(const MATRIX_T* W, MATRIX_T** Q, MATRIX_T** F) {
	if (safe_mode) check_status(W, SQUARE | WEIGHTED, "motif4struct_wei");
	
	// load motif34lib M4 M4n ID4 N4
	VECTOR_T* ID4;
	VECTOR_T* N4;
	MATRIX_T* M4 = motif4generate(&ID4, &N4);
	
	// n=length(W);
	int n = length(W);
	
	// I=zeros(199,n);
	MATRIX_T* I = zeros(199, n);
	
	// Q=zeros(199,n);
	if (Q != NULL) {
		*Q = zeros(199, n);
	}
	
	// F=zeros(199,n);
	if (F != NULL) {
		*F = zeros(199, n);
	}
	
	// A=1*(W~=0);
	MATRIX_T* A = compare_elements(W, fp_not_equal, 0.0);
	
	// As=A|A.';
	MATRIX_T* A_transpose = MATRIX_ID(alloc)(A->size2, A->size1);
	MATRIX_ID(transpose_memcpy)(A_transpose, A);
	MATRIX_T* As = logical_or(A, A_transpose);
	MATRIX_ID(free)(A_transpose);
	
	// for u=1:n-3
	for (int u = 0; u < n - 3; u++) {
		
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
				
				// V2=V2|([false(1,v1) As(u,v1+1:n)]);
				VECTOR_T* V2_1 = V2;
				VECTOR_T* V2_2 = VECTOR_ID(alloc)(n);
				MATRIX_ID(get_row)(V2_2, As, u);
				for (int i = 0; i <= v1; i++) {
					VECTOR_ID(set)(V2_2, i, 0.0);
				}
				V2 = logical_or(V2_1, V2_2);
				VECTOR_ID(free)(V2_1);
				VECTOR_ID(free)(V2_2);
				
				// for v2=find(V2)
				VECTOR_T* find_V2 = find(V2);
				if (find_V2 != NULL) {
					for (int i_find_V2 = 0; i_find_V2 < (int)find_V2->size; i_find_V2++) {
						int v2 = (int)VECTOR_ID(get)(find_V2, i_find_V2);
						
						// vz=max(v1,v2);
						int vz = (v1 > v2) ? v1 : v2;
						
						// V3=([false(1,u) As(v2,u+1:n)]);
						VECTOR_T* V3 = VECTOR_ID(alloc)(n);
						MATRIX_ID(get_row)(V3, As, v2);
						for (int i = 0; i <= u; i++) {
							VECTOR_ID(set)(V3, i, 0.0);
						}
						
						// V3(V2)=0;
						logical_index_assign(V3, V2, 0.0);
						
						// V3=V3|([false(1,v2) As(v1,v2+1:n)]);
						VECTOR_T* V3_1 = V3;
						VECTOR_T* V3_2 = VECTOR_ID(alloc)(n);
						MATRIX_ID(get_row)(V3_2, As, v1);
						for (int i = 0; i <= v2; i++) {
							VECTOR_ID(set)(V3_2, i, 0.0);
						}
						V3 = logical_or(V3_1, V3_2);
						VECTOR_ID(free)(V3_1);
						VECTOR_ID(free)(V3_2);
						
						// V3(V1)=0;
						logical_index_assign(V3, V1, 0.0);
						
						// V3=V3|([false(1,vz) As(u,vz+1:n)]);
						V3_1 = V3;
						V3_2 = VECTOR_ID(alloc)(n);
						MATRIX_ID(get_row)(V3_2, As, u);
						for (int i = 0; i <= vz; i++) {
							VECTOR_ID(set)(V3_2, i, 0.0);
						}
						V3 = logical_or(V3_1, V3_2);
						VECTOR_ID(free)(V3_1);
						VECTOR_ID(free)(V3_2);
						
						// for v3=find(V3)
						VECTOR_T* find_V3 = find(V3);
						if (find_V3 != NULL ) {
							for (int i_find_V3 = 0; i_find_V3 < (int)find_V3->size; i_find_V3++) {
								int v3 = (int)VECTOR_ID(get)(find_V3, i_find_V3);
								
								// w=[W(v1,u) W(v2,u) W(v3,u) W(u,v1) W(v2,v1) W(v3,v1) W(u,v2) W(v1,v2) W(v3,v2) W(u,v3) W(v1,v3) W(v2,v3)];
								int WA_rows[] = { v1, v2, v3, u, v2, v3, u, v1, v3, u, v1, v2 };
								int WA_cols[] = { u, u, u, v1, v1, v1, v2, v2, v2, v3, v3, v3 };
								VECTOR_T* w = VECTOR_ID(alloc)(12);
								for (int i = 0; i < 12; i++) {
									VECTOR_ID(set)(w, i, MATRIX_ID(get)(W, WA_rows[i], WA_cols[i]));
								}
								
								// s=uint64(sum(10.^(11:-1:0).*[A(v1,u) A(v2,u) A(v3,u) A(u,v1) A(v2,v1) A(v3,v1) A(u,v2) A(v1,v2) A(v3,v2) A(u,v3) A(v1,v3) A(v2,v3)]));
								VECTOR_T* s = VECTOR_ID(alloc)(12);
								for (int i = 0; i < 12; i++) {
									VECTOR_ID(set)(s, i, MATRIX_ID(get)(A, WA_rows[i], WA_cols[i]));
								}
								
								// ind=(s==M4n);
								int ind = 0;
								for ( ; ind < (int)M4->size1; ind++) {
									VECTOR_ID(view) M4_row_ind = MATRIX_ID(row)(M4, ind);
									if (compare_vectors(s, &M4_row_ind.vector) == 0) {
										break;
									}
								}
								VECTOR_ID(free)(s);
								if (ind < (int)M4->size1) {
									
									// M=w.*M4(ind,:);
									VECTOR_T* M = VECTOR_ID(alloc)(M4->size2);
									MATRIX_ID(get_row)(M, M4, ind);
									VECTOR_ID(mul)(M, w);
									
									// id=ID4(ind);
									int id = (int)VECTOR_ID(get)(ID4, ind) - 1;
									
									// l=N4(ind);
									int l = (int)VECTOR_ID(get)(N4, ind);
									
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
									
									// I(id,[u v1 v2 v3])=I(id,[u v1 v2 v3])+[i i i i];
									// Q(id,[u v1 v2 v3])=Q(id,[u v1 v2 v3])+[q q q q];
									// F(id,[u v1 v2 v3])=F(id,[u v1 v2 v3])+[1 1 1 1];
									int IQF_cols[] = { u, v1, v2, v3 };
									for (int j = 0; j < 4; j++) {
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
							
							VECTOR_ID(free)(find_V3);
						}
						
						VECTOR_ID(free)(V3);
					}
					
					VECTOR_ID(free)(find_V2);
				}
				
				VECTOR_ID(free)(V2);
			}
			
			VECTOR_ID(free)(find_V1);
		}
		
		VECTOR_ID(free)(V1);
	}
	
	VECTOR_ID(free)(ID4);
	VECTOR_ID(free)(N4);
	MATRIX_ID(free)(M4);
	MATRIX_ID(free)(A);
	MATRIX_ID(free)(As);
	return I;
}

/*
 * Returns per-motif metrics instead of per-motif, per-node metrics.
 */
VECTOR_T* BCT_NAMESPACE::motif4struct_wei_v(const MATRIX_T* W, VECTOR_T** Q, VECTOR_T** F) {
	MATRIX_T* _Q;
	MATRIX_T* _F;
	MATRIX_T* _I = motif4struct_wei(W, &_Q, &_F);
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
