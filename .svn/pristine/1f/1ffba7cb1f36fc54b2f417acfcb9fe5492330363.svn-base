#include "bct.h"

/*
 * Counts occurrences of four-node functional motifs in a binary graph.
 */
VECTOR_T* BCT_NAMESPACE::motif4funct_bin(const MATRIX_T* W, MATRIX_T** F) {
	if (safe_mode) check_status(W, SQUARE | BINARY, "motif4funct_bin");
	
	// load motif34lib M4 ID4 N4
	VECTOR_T* ID4;
	VECTOR_T* N4;
	MATRIX_T* M4 = motif4generate(&ID4, &N4);
	
	// n=length(W);
	int n = length(W);
	
	// f=zeros(199,1);
	VECTOR_T* f = zeros_vector(199);
	
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
								
								// a=[A(v1,u);A(v2,u);A(v3,u);A(u,v1);A(v2,v1);A(v3,v1);A(u,v2);A(v1,v2);A(v3,v2);A(u,v3);A(v1,v3);A(v2,v3)];
								int A_rows[] = { v1, v2, v3, u, v2, v3, u, v1, v3, u, v1, v2 };
								int A_cols[] = { u, u, u, v1, v1, v1, v2, v2, v2, v3, v3, v3 };
								MATRIX_T* a = MATRIX_ID(alloc)(12, 1);
								for (int i = 0; i < 12; i++) {
									MATRIX_ID(set)(a, i, 0, MATRIX_ID(get)(A, A_rows[i], A_cols[i]));
								}
								
								// ind=(M4*a)==N4;
								MATRIX_T* M4_mul_a_m = mul(M4, a);
								MATRIX_ID(free)(a);
								VECTOR_T* M4_mul_a = to_vector(M4_mul_a_m);
								MATRIX_ID(free)(M4_mul_a_m);
								VECTOR_T* ind = compare_elements(M4_mul_a, fp_equal, N4);
								VECTOR_ID(free)(M4_mul_a);
								
								// id=ID4(ind);
								VECTOR_T* id = logical_index(ID4, ind);
								VECTOR_ID(free)(ind);
								if (id != NULL) {
									VECTOR_ID(add_constant)(id, -1.0);
								
									// [idu j]=unique(id);
									VECTOR_T* j;
									VECTOR_T* idu = unique(id, "last", &j);
									VECTOR_ID(free)(id);
									
									// j=[0;j];
									VECTOR_T* temp = VECTOR_ID(alloc)(j->size + 1);
									VECTOR_ID(set)(temp, 0, -1.0);
									VECTOR_ID(view) temp_subv = VECTOR_ID(subvector)(temp, 1, j->size);
									VECTOR_ID(memcpy)(&temp_subv.vector, j);
									VECTOR_ID(free)(j);
									j = temp;
									
									// mu=length(idu);
									int mu = length(idu);
									
									// f2=zeros(mu,1);
									VECTOR_T* f2 = zeros_vector(mu);
									
									// for h=1:mu
									for (int h = 0; h < mu; h++) {
										
										// f2(h)=j(h+1)-j(h);
										FP_T j_h_add_1 = VECTOR_ID(get)(j, h + 1);
										FP_T j_h = VECTOR_ID(get)(j, h);
										VECTOR_ID(set)(f2, h, j_h_add_1 - j_h);
									}
									VECTOR_ID(free)(j);
									
									// f(idu)=f(idu)+f2;
									VECTOR_T* f_idu_add_f2 = ordinal_index(f, idu);
									VECTOR_ID(add)(f_idu_add_f2, f2);
									ordinal_index_assign(f, idu, f_idu_add_f2);
									VECTOR_ID(free)(f_idu_add_f2);
									
									// if nargout==2; F(idu,[u v1 v2 v3])=F(idu,[u v1 v2 v3])+[f2 f2 f2 f2]; end
									if (F != NULL) {
										FP_T F_cols[] = { (FP_T)u, (FP_T)v1, (FP_T)v2, (FP_T)v3 };
										VECTOR_ID(view) F_cols_vv = VECTOR_ID(view_array)(F_cols, 4);
										MATRIX_T* F_idx = ordinal_index(*F, idu, &F_cols_vv.vector);
										for (int i = 0; i < 4; i++) {
											VECTOR_ID(view) F_idx_col_i = MATRIX_ID(column)(F_idx, i);
											VECTOR_ID(add)(&F_idx_col_i.vector, f2);
										}
										ordinal_index_assign(*F, idu, &F_cols_vv.vector, F_idx);
										MATRIX_ID(free)(F_idx);
									}
									
									VECTOR_ID(free)(idu);
									VECTOR_ID(free)(f2);
								}
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
	return f;
}
