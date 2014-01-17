#include <cmath>

#include "bct.h"

/*
 * Finds paths from a set of source nodes up to a given length.  Note that there
 * is no savepths argument; if all paths are desired, pass a valid pointer as
 * the allpths argument.  There is also no tpath argument as its value may
 * overflow a C++ long.  Since 0 is a valid node index in C++, -1 is used as the
 * "filler" value in allpths rather than 0 as in MATLAB.  Pq (the main return),
 * plq, and util are indexed by path length.  They therefore have (qmax + 1)
 * elements and contain no valid data at index 0.
 */
std::vector<MATRIX_T*> BCT_NAMESPACE::findpaths(const MATRIX_T* CIJ, const VECTOR_T* sources, int qmax, VECTOR_T** plq, int* qstop, MATRIX_T** allpths, MATRIX_T** util) {
	if (safe_mode) check_status(CIJ, SQUARE, "findpaths");
	
	// CIJ = double(CIJ~=0);
	MATRIX_T* _CIJ = compare_elements(CIJ, fp_not_equal, 0.0);
	
	// N = size(CIJ,1);
	int N = CIJ->size1;
	
	// pths = [];
	MATRIX_T* pths = NULL;
	
	// Pq = zeros(N,N,qmax);
	std::vector<MATRIX_T*> Pq(qmax + 1);
	Pq[0] = NULL;
	for (int i = 1; i <= qmax; i++) {
		Pq[i] = zeros(N, N);
	}
	
	// util = zeros(N,qmax);
	if (util != NULL) {
		*util = zeros(N, qmax + 1);
	}
	
	// q = 1;
	int q = 1;
	
	VECTOR_T* _CIJ_cols = sequence(0, N - 1);
	MATRIX_T* _CIJ_idx = ordinal_index(_CIJ, sources, _CIJ_cols);
	VECTOR_ID(free)(_CIJ_cols);
	MATRIX_T* _CIJ_idx_ij = find_ij(_CIJ_idx);
	MATRIX_ID(free)(_CIJ_idx);
	pths = MATRIX_ID(alloc)(2, _CIJ_idx_ij->size1);
	for (int i = 0; i < (int)_CIJ_idx_ij->size1; i++) {
		int i_row = (int)MATRIX_ID(get)(_CIJ_idx_ij, i, 0);
		int i_start = (int)VECTOR_ID(get)(sources, i_row);
		int i_end = (int)MATRIX_ID(get)(_CIJ_idx_ij, i, 1);
		MATRIX_ID(set)(pths, 0, i, (FP_T)i_start);
		MATRIX_ID(set)(pths, 1, i, (FP_T)i_end);
	}
	MATRIX_ID(free)(_CIJ_idx_ij);
	
	// util(1:N,q) = util(1:N,q)+hist(reshape(pths,1,size(pths,1)*size(pths,2)),1:N)';
	if (util != NULL) {
		VECTOR_T* reshape_pths = to_vector(pths);
		VECTOR_T* centers = sequence(0, N - 1);
		VECTOR_T* hist_reshape_pths = hist(reshape_pths, centers);
		VECTOR_ID(free)(reshape_pths);
		VECTOR_ID(free)(centers);
		VECTOR_ID(view) util_col_q = MATRIX_ID(column)(*util, q);
		VECTOR_ID(add)(&util_col_q.vector, hist_reshape_pths);
		VECTOR_ID(free)(hist_reshape_pths);
	}
	
	// for np=1:size(pths,2)
	for (int np = 0; np < (int)pths->size2; np++) {
		
		// Pq(pths(1,np),pths(q+1,np),q) = Pq(pths(1,np),pths(q+1,np),q) + 1;
		int i = (int)MATRIX_ID(get)(pths, 0, np);
		int j = (int)MATRIX_ID(get)(pths, q, np);
		MATRIX_ID(set)(Pq[q], i, j, MATRIX_ID(get)(Pq[q], i, j) + 1.0);
	}
	
	// if (savepths==1)
	if (allpths != NULL) {
		
		// allpths = pths;
		*allpths = copy(pths);
	}
	
	// for q=2:qmax
	for (q = 2; q <= qmax; q++) {
		
		// npths = zeros(q+1,len_npths);
		// Handle as a std::vector<VECTOR_T*> rather than preallocating a MATRIX_T*
		std::vector<VECTOR_T*> npths_v;
		
		// endp = unique(pths(q,:));
		VECTOR_ID(view) pths_row_q_sub_1 = MATRIX_ID(row)(pths, q - 1);
		VECTOR_T* endp = unique(&pths_row_q_sub_1.vector);
		
		// npthscnt = 0;
		int npthscnt = 0;
		
		// for ii=1:length(endp)
		for (int ii = 0; ii < length(endp); ii++) {
			
			// i = endp(ii);
			int i = (int)VECTOR_ID(get)(endp, ii);
			
			// [pa,pb] = find(pths(q,:) == i);
			VECTOR_T* pths_row_q_sub_1_eq_i = compare_elements(&pths_row_q_sub_1.vector, fp_equal, (FP_T)i);
			VECTOR_T* pb = find(pths_row_q_sub_1_eq_i);
			VECTOR_ID(free)(pths_row_q_sub_1_eq_i);
			if (pb != NULL) {
				
				// nendp = find(CIJ(i,:)==1);
				VECTOR_ID(const_view) _CIJ_row_i = MATRIX_ID(const_row)(_CIJ, i);
				VECTOR_T* _CIJ_row_i_eq_1 = compare_elements(&_CIJ_row_i.vector, fp_equal, 1.0);
				VECTOR_T* nendp = find(_CIJ_row_i_eq_1);
				VECTOR_ID(free)(_CIJ_row_i_eq_1);
				
				// if (~isempty(nendp))
				if (nendp != NULL) {
					
					// for jj=1:length(nendp)
					for (int jj = 0; jj < length(nendp); jj++) {
						
						// j = nendp(jj);
						int j = (int)VECTOR_ID(get)(nendp, jj);
						
						// pb_temp = pb(sum(j==pths(2:q,pb),1)==0);
						VECTOR_T* pths_rows = sequence(1, q - 1);
						MATRIX_T* pths_idx = ordinal_index(pths, pths_rows, pb);
						VECTOR_ID(free)(pths_rows);
						MATRIX_T* pths_idx_eq_j = compare_elements(pths_idx, fp_equal, (FP_T)j);
						MATRIX_ID(free)(pths_idx);
						VECTOR_T* sum_pths_idx_eq_j = sum(pths_idx_eq_j, 1);
						MATRIX_ID(free)(pths_idx_eq_j);
						VECTOR_T* sum_pths_idx_eq_j_eq_0 = compare_elements(sum_pths_idx_eq_j, fp_equal, 0.0);
						VECTOR_ID(free)(sum_pths_idx_eq_j);
						VECTOR_T* pb_temp = logical_index(pb, sum_pths_idx_eq_j_eq_0);
						VECTOR_ID(free)(sum_pths_idx_eq_j_eq_0);
						if (pb_temp != NULL) {
						
							// npths(:,npthscnt+1:npthscnt+length(pb_temp)) = [pths(:,pb_temp)' ones(length(pb_temp),1)*j]';
							pths_rows = sequence(0, pths->size1 - 1);
							pths_idx = ordinal_index(pths, pths_rows, pb_temp);
							VECTOR_ID(free)(pths_rows);
							MATRIX_T* temp = MATRIX_ID(alloc)(pths_idx->size1 + 1, pths_idx->size2);
							MATRIX_ID(view) temp_subm = MATRIX_ID(submatrix)(temp, 0, 0, pths_idx->size1, pths_idx->size2);
							MATRIX_ID(memcpy)(&temp_subm.matrix, pths_idx);
							MATRIX_ID(free)(pths_idx);
							pths_idx = temp;
							VECTOR_T* last_row = VECTOR_ID(alloc)(length(pb_temp));
							VECTOR_ID(set_all)(last_row, (FP_T)j);
							MATRIX_ID(set_row)(pths_idx, pths_idx->size1 - 1, last_row);
							VECTOR_ID(free)(last_row);
							for (int i = 0; i < length(pb_temp); i++) {
								npths_v.push_back(zeros_vector(q + 1));
								VECTOR_ID(view) pths_idx_col_i = MATRIX_ID(column)(pths_idx, i);
								VECTOR_ID(view) npths_v_i_subv = VECTOR_ID(subvector)(npths_v[npthscnt + i], 0, pths_idx->size1);
								VECTOR_ID(memcpy)(&npths_v_i_subv.vector, &pths_idx_col_i.vector);
							}
							MATRIX_ID(free)(pths_idx);
							
							// npthscnt = npthscnt+length(pb_temp);
							npthscnt += length(pb_temp);
							
							// Pq(1:N,j,q) = Pq(1:N,j,q)+(hist(pths(1,pb_temp),1:N))';
							VECTOR_ID(view) pths_row_0 = MATRIX_ID(row)(pths, 0);
							VECTOR_T* pths_row_0_idx = ordinal_index(&pths_row_0.vector, pb_temp);
							VECTOR_ID(free)(pb_temp);
							VECTOR_T* centers = sequence(0, N - 1);
							VECTOR_T* hist_pths_idx = hist(pths_row_0_idx, centers);
							VECTOR_ID(free)(pths_row_0_idx);
							VECTOR_ID(free)(centers);
							VECTOR_ID(view) Pq_q_col_j = MATRIX_ID(column)(Pq[q], j);
							VECTOR_ID(add)(&Pq_q_col_j.vector, hist_pths_idx);
							VECTOR_ID(free)(hist_pths_idx);
						}
					}
					
					VECTOR_ID(free)(nendp);
				}
				
				VECTOR_ID(free)(pb);
			}
		}
		
		VECTOR_ID(free)(endp);
		MATRIX_T* npths = MATRIX_ID(alloc)(q + 1, npthscnt);
		for (int i = 0; i < npthscnt; i++) {
			MATRIX_ID(set_col)(npths, i, npths_v[i]);
			VECTOR_ID(free)(npths_v[i]);
		}
		
		// if (savepths==1)
		if (allpths != NULL) {
			
			// allpths = [allpths; zeros(1,size(allpths,2))];
			// allpths = [allpths npths(:,1:npthscnt)];
			MATRIX_T* temp = MATRIX_ID(alloc)((*allpths)->size1 + 1, (*allpths)->size2 + npthscnt);
			MATRIX_ID(set_all)(temp, -1.0);
			MATRIX_ID(view) temp_subm = MATRIX_ID(submatrix)(temp, 0, 0, (*allpths)->size1, (*allpths)->size2);
			MATRIX_ID(memcpy)(&temp_subm.matrix, *allpths);
			temp_subm = MATRIX_ID(submatrix)(temp, 0, (*allpths)->size2, (*allpths)->size1 + 1, npthscnt);
			MATRIX_ID(view) npths_subm = MATRIX_ID(submatrix)(npths, 0, 0, npths->size1, npthscnt);
			MATRIX_ID(memcpy)(&temp_subm.matrix, &npths_subm.matrix);
			MATRIX_ID(free)(*allpths);
			*allpths = temp;
		}
		
		// util(1:N,q) = util(1:N,q) + hist(reshape(npths,1,size(npths,1)*size(npths,2)),1:N)' - diag(Pq(:,:,q));
		if (util != NULL) {
			VECTOR_T* reshape_npths = to_vector(npths);
			VECTOR_T* centers = sequence(0, N - 1);
			VECTOR_T* hist_reshape_npths = hist(reshape_npths, centers);
			VECTOR_ID(free)(reshape_npths);
			VECTOR_ID(free)(centers);
			VECTOR_ID(view) util_col_q = MATRIX_ID(column)(*util, q);
			VECTOR_ID(add)(&util_col_q.vector, hist_reshape_npths);
			VECTOR_ID(free)(hist_reshape_npths);
			VECTOR_ID(view) diag_Pq_q = MATRIX_ID(diagonal)(Pq[q]);
			VECTOR_ID(sub)(&util_col_q.vector, &diag_Pq_q.vector);
		}
		
		// pths = npths(:,npths(1,:)~=npths(q+1,:));
		VECTOR_T* npths_rows = sequence(0, npths->size1 - 1);
		VECTOR_ID(view) npths_row_0 = MATRIX_ID(row)(npths, 0);
		VECTOR_ID(view) npths_row_q = MATRIX_ID(row)(npths, q);
		VECTOR_T* npths_cols = compare_elements(&npths_row_0.vector, fp_not_equal, &npths_row_q.vector);
		MATRIX_ID(free)(pths);
		pths = ord_log_index(npths, npths_rows, npths_cols);
		MATRIX_ID(free)(npths);
		VECTOR_ID(free)(npths_rows);
		VECTOR_ID(free)(npths_cols);
		
		// if (isempty(pths))
		// ...
		if (pths == NULL) {
			break;
		}
	}
	
	MATRIX_ID(free)(_CIJ);
	if (pths != NULL) {
		MATRIX_ID(free)(pths);
	}
	
	// qstop = q;
	if (qstop != NULL) {
		*qstop = (q <= qmax) ? q : qmax;
	}
	
	// tpath = sum(sum(sum(Pq)));
	// plq = reshape(sum(sum(Pq)),1,qmax)
	if (plq != NULL) {
		*plq = VECTOR_ID(alloc)(qmax + 1);
		VECTOR_ID(set)(*plq, 0, 0.0);
		for (int i = 1; i <= qmax; i++) {
			VECTOR_T* sum_Pq_i = sum(Pq[i]);
			VECTOR_ID(set)(*plq, i, sum(sum_Pq_i));
			VECTOR_ID(free)(sum_Pq_i);
		}
	}
	
	return Pq;
}
