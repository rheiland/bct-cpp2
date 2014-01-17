#include "bct.h"

/*
 * Computes the fraction of all paths that are cycles.
 */
VECTOR_T* BCT_NAMESPACE::cycprob_fcyc(const std::vector<MATRIX_T*>& Pq) {
	
	// fcyc = zeros(1,size(Pq,3));
	VECTOR_T* fcyc = zeros_vector(Pq.size());
	
	// for q=1:size(Pq,3)
	for (int q = 1; q < (int)Pq.size(); q++) {
		
		// if(sum(sum(Pq(:,:,q)))>0)
		VECTOR_T* sum_Pq_q = sum(Pq[q]);
		FP_T sum_sum_Pq_q = sum(sum_Pq_q);
		VECTOR_ID(free)(sum_Pq_q);
		if (sum_sum_Pq_q > 0.0) {
			
			// fcyc(q) = sum(diag(Pq(:,:,q)))/sum(sum(Pq(:,:,q)));
			VECTOR_ID(view) diag_Pq_q = MATRIX_ID(diagonal)(Pq[q]);
			FP_T sum_diag_Pq_q = sum(&diag_Pq_q.vector);
			VECTOR_ID(set)(fcyc, q, sum_diag_Pq_q / sum_sum_Pq_q);
		}
	}
	
	return fcyc;
}

/*
 * Computes the probability that a non-cyclic path of length (q - 1) can be
 * extended to form a cycle of length q.
 */
VECTOR_T* BCT_NAMESPACE::cycprob_pcyc(const std::vector<MATRIX_T*>& Pq) {
	
	// pcyc = zeros(1,size(Pq,3));
	VECTOR_T* pcyc = zeros_vector(Pq.size());
	
	// for q=2:size(Pq,3)
	for (int q = 2; q < (int)Pq.size(); q++) {
		
		// if((sum(sum(Pq(:,:,q-1)))-sum(diag(Pq(:,:,q-1))))>0)
		VECTOR_T* sum_Pq_q_sub_1 = sum(Pq[q - 1]);
		FP_T sum_sum_Pq_q_sub_1 = sum(sum_Pq_q_sub_1);
		VECTOR_ID(free)(sum_Pq_q_sub_1);
		VECTOR_ID(view) diag_Pq_q_sub_1 = MATRIX_ID(diagonal)(Pq[q - 1]);
		FP_T sum_diag_Pq_q_sub_1 = sum(&diag_Pq_q_sub_1.vector);
		if (sum_sum_Pq_q_sub_1 - sum_diag_Pq_q_sub_1 > 0.0) {
			
			// pcyc(q) = sum(diag(Pq(:,:,q)))/(sum(sum(Pq(:,:,q-1)))-sum(diag(Pq(:,:,q-1))));
			VECTOR_ID(view) diag_Pq_q = MATRIX_ID(diagonal)(Pq[q]);
			FP_T sum_diag_Pq_q = sum(&diag_Pq_q.vector);
			VECTOR_ID(set)(pcyc, q, sum_diag_Pq_q / (sum_sum_Pq_q_sub_1 - sum_diag_Pq_q_sub_1));
		}
	}
	
	return pcyc;
}
