#include "bct.h"

/*
 * Generates a random directed binary graph with a ring lattice organization.
 */
MATRIX_T* BCT_NAMESPACE::makeringlatticeCIJ(int N, int K) {
	
	// CIJ = zeros(N);
	MATRIX_T* CIJ = zeros(N);
	
	// CIJ1 = ones(N);
	MATRIX_T* CIJ1 = ones(N);
	
	// KK = 0;
	int KK = 0;
	
	// cnt = 0;
	int cnt = -1;
	
	// seq = 1:N-1;
	VECTOR_T* seq = sequence(1, N - 1);
	
	// seq2 = N-1:-1:1;
	VECTOR_T* seq2 = sequence(N - 1, -1, 1);
	
	MATRIX_T* dCIJ = NULL;
	
	// while (KK<K)
	while (KK < K) {
		
		// cnt = cnt + 1;
		cnt++;
		
		// dCIJ = triu(CIJ1,seq(cnt))-triu(CIJ1,seq(cnt)+1);
		if (dCIJ != NULL) {
			MATRIX_ID(free)(dCIJ);
		}
		dCIJ = triu(CIJ1, (int)VECTOR_ID(get)(seq, cnt));
		MATRIX_T* triu_CIJ1_seq_cnt_add_1 = triu(CIJ1, (int)VECTOR_ID(get)(seq, cnt) + 1);
		MATRIX_ID(sub)(dCIJ, triu_CIJ1_seq_cnt_add_1);
		MATRIX_ID(free)(triu_CIJ1_seq_cnt_add_1);
		
		// dCIJ2 = triu(CIJ1,seq2(cnt))-triu(CIJ1,seq2(cnt)+1);
		MATRIX_T* dCIJ2 = triu(CIJ1, (int)VECTOR_ID(get)(seq2, cnt));
		MATRIX_T* triu_CIJ1_seq2_cnt_add_1 = triu(CIJ1, (int)VECTOR_ID(get)(seq2, cnt) + 1);
		MATRIX_ID(sub)(dCIJ2, triu_CIJ1_seq2_cnt_add_1);
		MATRIX_ID(free)(triu_CIJ1_seq2_cnt_add_1);
		
		// dCIJ = dCIJ+dCIJ'+dCIJ2+dCIJ2';
		MATRIX_T* dCIJ_transpose = MATRIX_ID(alloc)(dCIJ->size2, dCIJ->size1);
		MATRIX_ID(transpose_memcpy)(dCIJ_transpose, dCIJ);
		MATRIX_T* dCIJ2_transpose = MATRIX_ID(alloc)(dCIJ2->size2, dCIJ2->size1);
		MATRIX_ID(transpose_memcpy)(dCIJ2_transpose, dCIJ2);
		MATRIX_ID(add)(dCIJ, dCIJ_transpose);
		MATRIX_ID(free)(dCIJ_transpose);
		MATRIX_ID(add)(dCIJ, dCIJ2);
		MATRIX_ID(free)(dCIJ2);
		MATRIX_ID(add)(dCIJ, dCIJ2_transpose);
		MATRIX_ID(free)(dCIJ2_transpose);
		
		// CIJ = CIJ + dCIJ;
		MATRIX_ID(add)(CIJ, dCIJ);
		
		// KK = sum(sum(CIJ));
		VECTOR_T* sum_CIJ = sum(CIJ);
		KK = (int)sum(sum_CIJ);
		VECTOR_ID(free)(sum_CIJ);
	}
	
	MATRIX_ID(free)(CIJ1);
	VECTOR_ID(free)(seq);
	VECTOR_ID(free)(seq2);
	
	// overby = KK-K;
	int overby = KK - K;
	
	// if(overby>0)
	if (overby > 0) {
		
		// [i j] = find(dCIJ);
		MATRIX_T* find_dCIJ = find_ij(dCIJ);
		VECTOR_ID(view) i = MATRIX_ID(column)(find_dCIJ, 0);
		VECTOR_ID(view) j = MATRIX_ID(column)(find_dCIJ, 1);
		
		// rp = randperm(length(i));
		gsl_permutation* rp = randperm(length(&i.vector));
		
		// for ii=1:overby
		for (int ii = 0; ii < overby; ii++) {
			
			// CIJ(i(rp(ii)),j(rp(ii))) = 0;
			int i_rp_ii = VECTOR_ID(get)(&i.vector, gsl_permutation_get(rp, ii));
			int j_rp_ii = VECTOR_ID(get)(&j.vector, gsl_permutation_get(rp, ii));
			MATRIX_ID(set)(CIJ, i_rp_ii, j_rp_ii, 0.0);
		}
		
		MATRIX_ID(free)(find_dCIJ);
		gsl_permutation_free(rp);
	}
	
	MATRIX_ID(free)(dCIJ);
	return CIJ;
}
