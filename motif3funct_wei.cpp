#include "bct.h"

/*
 * Counts occurrences of three-node functional motifs in a weighted graph.
 * Returns intensity and (optionally) coherence and motif counts.
 */
MATRIX_T* BCT_NAMESPACE::motif3funct_wei(const MATRIX_T* W, MATRIX_T** Q, MATRIX_T** F) {
	if (safe_mode) check_status(W, SQUARE | WEIGHTED, "motif3funct_wei");
	
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
						
						// a=[A(v1,u);A(v2,u);A(u,v1);A(v2,v1);A(u,v2);A(v1,v2)];
						MATRIX_T* a = MATRIX_ID(alloc)(6, 1);
						for (int i = 0; i < 6; i++) {
							MATRIX_ID(set)(a, i, 0, MATRIX_ID(get)(A, WA_rows[i], WA_cols[i]));
						}
						
						// ind=(M3*a)==N3;
						MATRIX_T* M3_mul_a_m = mul(M3, a);
						MATRIX_ID(free)(a);
						VECTOR_T* M3_mul_a = to_vector(M3_mul_a_m);
						MATRIX_ID(free)(M3_mul_a_m);
						VECTOR_T* ind = compare_elements(M3_mul_a, fp_equal, N3);
						VECTOR_ID(free)(M3_mul_a);
						
						// m=sum(ind);
						int m = (int)sum(ind);
						if (m > 0) {
							
							// M=M3(ind,:).*repmat(w,m,1);
							VECTOR_T* M3_cols = sequence(0, M3->size2 - 1);
							MATRIX_T* M = log_ord_index(M3, ind, M3_cols);
							VECTOR_ID(free)(M3_cols);
							for (int i = 0; i < (int)M->size1; i++) {
								VECTOR_ID(view) M_row_i = MATRIX_ID(row)(M, i);
								VECTOR_ID(mul)(&M_row_i.vector, w);
							}
							
							// id=ID3(ind);
							VECTOR_T* id = logical_index(ID3, ind);
							VECTOR_ID(add_constant)(id, -1.0);
							
							// l=N3(ind);
							VECTOR_T* l = logical_index(N3, ind);
							
							// x=sum(M,2)./l;
							VECTOR_T* x = sum(M, 2);
							VECTOR_ID(div)(x, l);
							
							// M(M==0)=1;
							MATRIX_T* M_eq_0 = compare_elements(M, fp_equal, 0.0);
							logical_index_assign(M, M_eq_0, 1.0);
							MATRIX_ID(free)(M_eq_0);
							
							// i=prod(M,2).^(1./l);
							VECTOR_T* prod_M = prod(M, 2);
							MATRIX_ID(free)(M);
							VECTOR_T* l_pow_neg_1 = pow_elements(l, -1.0);
							VECTOR_ID(free)(l);
							VECTOR_T* i = pow_elements(prod_M, l_pow_neg_1);
							VECTOR_ID(free)(prod_M);
							VECTOR_ID(free)(l_pow_neg_1);
							
							// q = i./x;
							VECTOR_T* q = copy(i);
							VECTOR_ID(div)(q, x);
							VECTOR_ID(free)(x);
							
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
							
							// i2=zeros(mu,1);
							VECTOR_T* i2 = zeros_vector(mu);
							
							// q2=i2; f2=i2;
							VECTOR_T* q2 = copy(i2);
							VECTOR_T* f2 = copy(i2);
							
							// for h=1:mu
							for (int h = 0; h < mu; h++) {
								
								// i2(h)=sum(i(j(h)+1:j(h+1)));
								int j_h = (int)VECTOR_ID(get)(j, h);
								int j_h_add_1 = (int)VECTOR_ID(get)(j, h + 1);
								VECTOR_T* iq_indices = sequence(j_h + 1, j_h_add_1);
								VECTOR_T* i_idx = ordinal_index(i, iq_indices);
								VECTOR_ID(set)(i2, h, sum(i_idx));
								VECTOR_ID(free)(i_idx);
								
								// q2(h)=sum(q(j(h)+1:j(h+1)));
								VECTOR_T* q_idx = ordinal_index(q, iq_indices);
								VECTOR_ID(free)(iq_indices);
								VECTOR_ID(set)(q2, h, sum(q_idx));
								VECTOR_ID(free)(q_idx);
								
								// f2(h)=j(h+1)-j(h);
								VECTOR_ID(set)(f2, h, (FP_T)(j_h_add_1 - j_h));
							}
							VECTOR_ID(free)(i);
							VECTOR_ID(free)(q);
							VECTOR_ID(free)(j);
							
							// I(idu,[u v1 v2])=I(idu,[u v1 v2])+[i2 i2 i2];
							// Q(idu,[u v1 v2])=Q(idu,[u v1 v2])+[q2 q2 q2];
							// F(idu,[u v1 v2])=F(idu,[u v1 v2])+[f2 f2 f2];
							FP_T IQF_cols[] = { (FP_T)u, (FP_T)v1, (FP_T)v2 };
							VECTOR_ID(view) IQF_cols_vv = VECTOR_ID(view_array)(IQF_cols, 3);
							MATRIX_T* I_idx = ordinal_index(I, idu, &IQF_cols_vv.vector);
							MATRIX_T* Q_idx = NULL;
							MATRIX_T* F_idx = NULL;
							if (Q != NULL) {
								Q_idx = ordinal_index(*Q, idu, &IQF_cols_vv.vector);
							}
							if (F != NULL) {
								F_idx = ordinal_index(*F, idu, &IQF_cols_vv.vector);
							}
							for (int j = 0; j < 3; j++) {
								VECTOR_ID(view) I_idx_col_j = MATRIX_ID(column)(I_idx, j);
								VECTOR_ID(add)(&I_idx_col_j.vector, i2);
								if (Q != NULL) {
									VECTOR_ID(view) Q_idx_col_j = MATRIX_ID(column)(Q_idx, j);
									VECTOR_ID(add)(&Q_idx_col_j.vector, q2);
								}
								if (F != NULL) {
									VECTOR_ID(view) F_idx_col_j = MATRIX_ID(column)(F_idx, j);
									VECTOR_ID(add)(&F_idx_col_j.vector, f2);
								}
							}
							VECTOR_ID(free)(i2);
							VECTOR_ID(free)(q2);
							VECTOR_ID(free)(f2);
							ordinal_index_assign(I, idu, &IQF_cols_vv.vector, I_idx);
							MATRIX_ID(free)(I_idx);
							if (Q != NULL) {
								ordinal_index_assign(*Q, idu, &IQF_cols_vv.vector, Q_idx);
								MATRIX_ID(free)(Q_idx);
							}
							if (F != NULL) {
								ordinal_index_assign(*F, idu, &IQF_cols_vv.vector, F_idx);
								MATRIX_ID(free)(F_idx);
							}
							VECTOR_ID(free)(idu);
						}
						
						VECTOR_ID(free)(w);
						VECTOR_ID(free)(ind);
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
VECTOR_T* BCT_NAMESPACE::motif3funct_wei_v(const MATRIX_T* W, VECTOR_T** Q, VECTOR_T** F) {
	MATRIX_T* _Q;
	MATRIX_T* _F;
	MATRIX_T* _I = motif3funct_wei(W, &_Q, &_F);
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
