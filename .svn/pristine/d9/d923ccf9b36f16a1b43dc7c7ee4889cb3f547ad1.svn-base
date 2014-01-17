#include "bct.h"

/*
 * Returns all motif isomorphs for a given motif ID and size.
 */
std::vector<MATRIX_T*> BCT_NAMESPACE::find_motif34(int m, int n) {
	
	// if n==3
	if (n == 3) {
		
		// load motif34lib M3 ID3
		VECTOR_T* ID3;
		MATRIX_T* M3 = motif3generate(&ID3);
		
		// ind=find(ID3==m).';
		VECTOR_T* ID3_eq_m = compare_elements(ID3, fp_equal, (FP_T)m);
		VECTOR_ID(free)(ID3);
		VECTOR_T* ind = find(ID3_eq_m);
		VECTOR_ID(free)(ID3_eq_m);
		
		// M=zeros(3,3,length(ind));
		std::vector<MATRIX_T*> M(length(ind));
		
		int i_nondiag[] = { 1, 2, 3, 5, 6, 7 };
		
		// for i=1:length(ind)
		for (int i = 0; i < length(ind); i++) {
			
			// M(:,:,i)=reshape([0 M3(ind(i),1:3) 0 M3(ind(i),4:6) 0],3,3);
			M[i] = MATRIX_ID(calloc)(3, 3);
			int ind_i = (int)VECTOR_ID(get)(ind, i);
			for (int j = 0; j < 6; j++) {
				ordinal_index_assign(M[i], i_nondiag[j], MATRIX_ID(get)(M3, ind_i, j));
			}
		}
		
		MATRIX_ID(free)(M3);
		VECTOR_ID(free)(ind);
		return M;
		
		// elseif n==4
	} else if (n == 4) {
		
		// load motif34lib M4 ID4;
		VECTOR_T* ID4;
		MATRIX_T* M4 = motif4generate(&ID4);
		
		// ind=find(ID4==m).';
		VECTOR_T* ID4_eq_m = compare_elements(ID4, fp_equal, (FP_T)m);
		VECTOR_ID(free)(ID4);
		VECTOR_T* ind = find(ID4_eq_m);
		VECTOR_ID(free)(ID4_eq_m);
		
		// M=zeros(4,4,length(ind));
		std::vector<MATRIX_T*> M(length(ind));
		
		int i_nondiag[] = { 1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14 };
		
		// for i=1:length(ind)
		for (int i = 0; i < length(ind); i++) {
			
			// M(:,:,i)=reshape([0 M4(ind(i),1:4) 0 M4(ind(i),5:8) 0 M4(ind(i),9:12) 0],4,4);
			M[i] = MATRIX_ID(calloc)(4, 4);
			int ind_i = (int)VECTOR_ID(get)(ind, i);
			for (int j = 0; j < 12; j++) {
				ordinal_index_assign(M[i], i_nondiag[j], MATRIX_ID(get)(M4, ind_i, j));
			}
		}
		
		MATRIX_ID(free)(M4);
		VECTOR_ID(free)(ind);
		return M;
	} else {
		return std::vector<MATRIX_T*>();
	}
}

/*
 * Returns the motif ID for a given matrix.
 */
int BCT_NAMESPACE::find_motif34(const MATRIX_T* m) {
	if (safe_mode) check_status(m, SQUARE, "find_motif34");
	
	// n=size(m,1);
	int n = m->size1;
	
	// M=eval(['find(motif' int2str(n) 'struct_bin(m))']);
	VECTOR_T* f;
	if (n == 3) {
		f = motif3struct_bin(m);
	} else if (n == 4) {
		f = motif4struct_bin(m);
	} else {
		return 0;
	}
	VECTOR_T* M_v = find(f);
	VECTOR_ID(free)(f);
	if (M_v != NULL) {
		int M = (int)VECTOR_ID(get)(M_v, 0);
		VECTOR_ID(free)(M_v);
		return M + 1;
	} else {
		return 0;
	}
}
