#include "bct.h"

/*
 * Our implementation of the BCT motif library does not include Mn ("M as a
 * single number") because a C++ long is generally not large enough to contain
 * the 12-digit numbers required for the four-node motif library.
 */

BCT_NAMESPACE::motif_mode_enum BCT_NAMESPACE::motif_mode = MILO;

BCT_NAMESPACE::motif_mode_enum BCT_NAMESPACE::get_motif_mode() { return motif_mode; }
void BCT_NAMESPACE::set_motif_mode(motif_mode_enum motif_mode) { BCT_NAMESPACE::motif_mode = motif_mode; }

/*
 * Constructs the three-node motif library.
 */
MATRIX_T* BCT_NAMESPACE::motif3generate(VECTOR_T** ID, VECTOR_T** N) {
	static MATRIX_T* M = NULL;
	static VECTOR_T* _ID = NULL;
	static VECTOR_T* _N = NULL;
	if (M == NULL) {
		
		// n=0;
		int n = -1;
		
		// M=false(54,6);
		M = MATRIX_ID(calloc)(54, 6);
		
		// CL=zeros(54,6,'uint8');
		MATRIX_T* CL = zeros(54, 6);
		
		// cl=zeros(1,6,'uint8');
		VECTOR_T* cl = zeros_vector(6);
		
		FP_T i_nondiag[] = { 1, 2, 3, 5, 6, 7 };
		VECTOR_ID(view) i_nondiag_vv = VECTOR_ID(view_array)(i_nondiag, 6);
		
		// for i=0:2^6-1
		for (int i = 0; i < 64; i++) {
			
			// m=dec2bin(i);
			// m=[num2str(zeros(1,6-length(m)), '%d') m];
			std::string m = dec2bin(i, 6);
			
			// G=str2num ([ ...
			// '0'   ' '  m(3)  ' '  m(5) ;
			// m(1)  ' '  '0'   ' '  m(6) ;
			// m(2)  ' '  m(4)  ' '  '0'   ]);
			MATRIX_T* G = MATRIX_ID(calloc)(3, 3);
			for (int i = 0; i < 6; i++) {
				int index = (int)i_nondiag[i];
				if (m[i] == '1') {
					ordinal_index_assign(G, index, 1.0);
				}
			}
			
			// Ko=sum(G,2);
			VECTOR_T* Ko = sum(G, 2);
			
			// Ki=sum(G,1).';
			VECTOR_T* Ki = sum(G, 1);
			
			// if Ko+Ki,
			VECTOR_T* Ko_add_Ki = copy(Ko);
			VECTOR_ID(add)(Ko_add_Ki, Ki);
			bool Ko_add_Ki_bool = to_bool(Ko_add_Ki);
			VECTOR_ID(free)(Ko_add_Ki);
			if (Ko_add_Ki_bool) {
				
				// n=n+1;
				n++;
				
				// cl(:)=sortrows([Ko Ki]).';
				MATRIX_T* Ko_Ki = concatenate_rows(Ko, Ki);
				MATRIX_T* Ko_Ki_sorted = sortrows(Ko_Ki);
				MATRIX_ID(free)(Ko_Ki);
				MATRIX_T* Ko_Ki_transpose = MATRIX_ID(alloc)(2, G->size1);
				MATRIX_ID(transpose_memcpy)(Ko_Ki_transpose, Ko_Ki_sorted);
				MATRIX_ID(free)(Ko_Ki_sorted);
				VECTOR_ID(free)(cl);
				cl = to_vector(Ko_Ki_transpose);
				MATRIX_ID(free)(Ko_Ki_transpose);
				
				// CL(n,:)=cl;
				MATRIX_ID(set_row)(CL, n, cl);
				
				// M(n,:)=G([2:4 6:8]);
				VECTOR_T* G_nondiag = ordinal_index(G, &i_nondiag_vv.vector);
				MATRIX_ID(set_row)(M, n, G_nondiag);
				VECTOR_ID(free)(G_nondiag);
			}
			
			MATRIX_ID(free)(G);
			VECTOR_ID(free)(Ko);
			VECTOR_ID(free)(Ki);
		}
		
		VECTOR_ID(free)(cl);
		
		// [u1 u2 ID]=unique(CL,'rows');
		MATRIX_T* u1 = unique_rows(CL, "last", NULL, &_ID);
		MATRIX_ID(free)(CL);
		MATRIX_ID(free)(u1);
		VECTOR_ID(add_constant)(_ID, 1.0);

		// id_mika=  [1  3  4  6  7  8  11];
		int id_mika[] = { 1, 3, 4, 6, 7, 8, 11 };
		
		// id_olaf= -[3  6  1 11  4  7   8];
		int id_olaf[] = { 3, 6, 1, 11, 4, 7, 8};
		
		// %convert IDs into Sporns & Kotter classification
		if (motif_mode == SPORNS) {
			for (int i = 0; i < (int)_ID->size; i++) {
				for (int j = 0; j < 7; j++) {
					if ((int)VECTOR_ID(get)(_ID, i) == id_mika[j]) {
						VECTOR_ID(set)(_ID, i, id_olaf[j]);
						break;
					}
				}
			}
		}
		
		// [X ind]=sortrows(ID);
		VECTOR_T* ind_v;
		VECTOR_T* X = sortrows(_ID, &ind_v);
		VECTOR_ID(free)(X);
		gsl_permutation* ind = to_permutation(ind_v);
		VECTOR_ID(free)(ind_v);
		
		// ID=ID(ind,:);
		VECTOR_T* _ID_permuted = permute(ind, _ID);
		VECTOR_ID(free)(_ID);
		_ID = _ID_permuted;
		
		// M=M(ind,:);
		MATRIX_T* M_permuted = permute_rows(ind, M);
		MATRIX_ID(free)(M);
		M = M_permuted;
		gsl_permutation_free(ind);
		
		// N=sum(M,2);
		_N = sum(M, 2);
	}
	
	if (ID != NULL) {
		*ID = copy(_ID);
	}
	if (N != NULL) {
		*N = copy(_N);
	}
	return copy(M);
}

/*
 * Constructs the four-node motif library.
 */
MATRIX_T* BCT_NAMESPACE::motif4generate(VECTOR_T** ID, VECTOR_T** N) {
	static MATRIX_T* M = NULL;
	static VECTOR_T* _ID = NULL;
	static VECTOR_T* _N = NULL;
	if (M == NULL) {
		
		// n=0;
		int n = -1;
		
		// M=false(3834,12);
		M = MATRIX_ID(calloc)(3834, 12);
		
		// CL=zeros(3834,16,'uint8');
		MATRIX_T* CL = zeros(3834, 16);
		
		// cl=zeros(1,16,'uint8');
		VECTOR_T* cl = zeros_vector(16);
		
		FP_T i_nondiag[] = { 1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14 };
		VECTOR_ID(view) i_nondiag_vv = VECTOR_ID(view_array)(i_nondiag, 12);
		
		// for i=0:2^12-1
		for (int i = 0; i < 4096; i++) {
			
			// m=dec2bin(i);
			// m=[num2str(zeros(1,12-length(m)), '%d') m];
			std::string m = dec2bin(i, 12);
			
			// G=str2num ([ ...
			// '0'   ' '  m(4)  ' '  m(7)  ' '  m(10) ;
			// m(1)  ' '  '0'   ' '  m(8)  ' '  m(11) ;
			// m(2)  ' '  m(5)  ' '  '0'   ' '  m(12) ;
			// m(3)  ' '  m(6)  ' '  m(9)  ' '  '0'    ]);
			MATRIX_T* G = MATRIX_ID(calloc)(4, 4);
			for (int i = 0; i < 12; i++) {
				int index = (int)i_nondiag[i];
				if (m[i] == '1') {
					ordinal_index_assign(G, index, 1.0);
				}
			}
			
			// Gs=G+G.';
			MATRIX_T* Gs = MATRIX_ID(alloc)(4, 4);
			MATRIX_ID(transpose_memcpy)(Gs, G);
			MATRIX_ID(add)(Gs, G);
			
			// v=Gs(1,:);
			VECTOR_T* v = VECTOR_ID(alloc)(4);
			MATRIX_ID(get_row)(v, Gs, 0);
			
			// for j=1:2,
			for (int j = 1; j <= 2; j++) {
				
				// v=any(Gs(v~=0,:),1)+v;
				VECTOR_T* v_neq_0 = compare_elements(v, fp_not_equal, 0.0);
				VECTOR_T* Gs_cols = sequence(0, Gs->size2 - 1);
				MATRIX_T* Gs_idx = log_ord_index(Gs, v_neq_0, Gs_cols);
				VECTOR_ID(free)(v_neq_0);
				VECTOR_ID(free)(Gs_cols);
				if (Gs_idx != NULL) {
					VECTOR_T* any_Gs_idx = any(Gs_idx, 1);
					VECTOR_ID(add)(v, any_Gs_idx);
					MATRIX_ID(free)(Gs_idx);
					VECTOR_ID(free)(any_Gs_idx);
				}
			}
			
			MATRIX_ID(free)(Gs);
			
			// if v
			bool v_bool = to_bool(v);
			VECTOR_ID(free)(v);
			if (v_bool) {
				
				// n=n+1;
				n++;
				
				// G2=(G*G)~=0;
				MATRIX_T* G_mul_G = mul(G, G);
				MATRIX_T* G2 = compare_elements(G_mul_G, fp_not_equal, 0.0);
				MATRIX_ID(free)(G_mul_G);
				
				// Ko=sum(G,2);
				VECTOR_T* Ko = sum(G, 2);
				
				// Ki=sum(G,1).';
				VECTOR_T* Ki = sum(G, 1);
				
				// Ko2=sum(G2,2);
				VECTOR_T* Ko2 = sum(G2, 2);
				
				// Ki2=sum(G2,1).';
				VECTOR_T* Ki2 = sum(G2, 1);
				
				// cl(:)=sortrows([Ki Ko Ki2 Ko2]).';
				MATRIX_T* Ki_Ko = concatenate_rows(Ki, Ko);
				VECTOR_ID(free)(Ki);
				VECTOR_ID(free)(Ko);
				MATRIX_T* Ki_Ko_Ki2 = concatenate_rows(Ki_Ko, Ki2);
				VECTOR_ID(free)(Ki2);
				MATRIX_ID(free)(Ki_Ko);
				MATRIX_T* Ki_Ko_Ki2_Ko2 = concatenate_rows(Ki_Ko_Ki2, Ko2);
				VECTOR_ID(free)(Ko2);
				MATRIX_ID(free)(Ki_Ko_Ki2);
				MATRIX_T* Ks_sorted = sortrows(Ki_Ko_Ki2_Ko2);
				MATRIX_ID(free)(Ki_Ko_Ki2_Ko2);
				MATRIX_T* Ks_transpose = MATRIX_ID(alloc)(4, G->size1);
				MATRIX_ID(transpose_memcpy)(Ks_transpose, Ks_sorted);
				MATRIX_ID(free)(Ks_sorted);
				VECTOR_ID(free)(cl);
				cl = to_vector(Ks_transpose);
				MATRIX_ID(free)(Ks_transpose);
				
				// CL(n,:)=cl;
				MATRIX_ID(set_row)(CL, n, cl);
				
				// M(n,:)=G([2:5 7:10 12:15]);
				VECTOR_T* G_nondiag = ordinal_index(G, &i_nondiag_vv.vector);
				MATRIX_ID(set_row)(M, n, G_nondiag);
				VECTOR_ID(free)(G_nondiag);
			}
			
			MATRIX_ID(free)(G);
		}
		
		VECTOR_ID(free)(cl);
		
		// [u1 u2 ID]=unique(CL,'rows');
		MATRIX_T* u1 = unique_rows(CL, "last", NULL, &_ID);
		MATRIX_ID(free)(CL);
		MATRIX_ID(free)(u1);
		VECTOR_ID(add_constant)(_ID, 1.0);
		
		// [X ind]=sortrows(ID);
		VECTOR_T* ind_v;
		VECTOR_T* X = sortrows(_ID, &ind_v);
		VECTOR_ID(free)(X);
		gsl_permutation* ind = to_permutation(ind_v);
		VECTOR_ID(free)(ind_v);
		
		// ID=ID(ind,:);
		VECTOR_T* _ID_permuted = permute(ind, _ID);
		VECTOR_ID(free)(_ID);
		_ID = _ID_permuted;
		
		// M=M(ind,:);
		MATRIX_T* M_permuted = permute_rows(ind, M);
		MATRIX_ID(free)(M);
		M = M_permuted;
		gsl_permutation_free(ind);
		
		// N=sum(M,2);
		_N = sum(M, 2);
	}
	
	if (ID != NULL) {
		*ID = copy(_ID);
	}
	if (N != NULL) {
		*N = copy(_N);
	}
	return copy(M);
}
