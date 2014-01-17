#include "bct.h"

/*
 * Counts occurrences of three-node functional motifs in a binary graph.
 */
VECTOR_T* BCT_NAMESPACE::motif3funct_bin(const MATRIX_T* W, MATRIX_T** F) {
	if (safe_mode) check_status(W, SQUARE | BINARY, "motif3funct_bin");
	
	// load motif34lib M3 ID3 N3
	VECTOR_T* ID3;
	VECTOR_T* N3;
	MATRIX_T* M3 = motif3generate(&ID3, &N3);
	
	// n=length(W);
	int n = length(W);
	
	// f=zeros(13,1);
	VECTOR_T* f = zeros_vector(13);
	
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
						
						// a=[A(v1,u);A(v2,u);A(u,v1);A(v2,v1);A(u,v2);A(v1,v2)];
						int A_rows[] = { v1, v2, u, v2, u, v1 };
						int A_cols[] = { u, u, v1, v1, v2, v2 };
						MATRIX_T* a = MATRIX_ID(alloc)(6, 1);
						for (int i = 0; i < 6; i++) {
							MATRIX_ID(set)(a, i, 0, MATRIX_ID(get)(A, A_rows[i], A_cols[i]));
						}
						
						// ind=(M3*a)==N3;
						MATRIX_T* M3_mul_a_m = mul(M3, a);
						MATRIX_ID(free)(a);
						VECTOR_T* M3_mul_a = to_vector(M3_mul_a_m);
						MATRIX_ID(free)(M3_mul_a_m);
						VECTOR_T* ind = compare_elements(M3_mul_a, fp_equal, N3);
						VECTOR_ID(free)(M3_mul_a);
						
						// id=ID3(ind);
						VECTOR_T* id = logical_index(ID3, ind);
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
							
							// if nargout==2; F(idu,[u v1 v2])=F(idu,[u v1 v2])+[f2 f2 f2]; end
							if (F != NULL) {
								FP_T F_cols[] = { (FP_T)u, (FP_T)v1, (FP_T)v2 };
								VECTOR_ID(view) F_cols_vv = VECTOR_ID(view_array)(F_cols, 3);
								MATRIX_T* F_idx = ordinal_index(*F, idu, &F_cols_vv.vector);
								for (int i = 0; i < 3; i++) {
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
	return f;
}
