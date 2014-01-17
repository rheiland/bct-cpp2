#include <cmath>

#include "bct.h"

/*
 * Generates a random directed binary graph with the given in-degree and out-
 * degree sequences.  Returns NULL if the algorithm failed to generate a graph
 * satisfying the given degree sequences.
 */
MATRIX_T* BCT_NAMESPACE::makerandCIJdegreesfixed(const VECTOR_T* in, const VECTOR_T* out) {
	gsl_rng* rng = get_rng();
	
	// n = length(in);
	int n = length(in);
	
	// k = sum(in);
	int k = sum(in);
	
	// inInv = zeros(k,1);
	VECTOR_T* inInv = zeros_vector(k);
	
	// outInv = inInv;
	VECTOR_T* outInv = copy(inInv);
	
	// iIn = 1; iOut = 1;
	int iIn = 0;
	int iOut = 0;
	
	// for i = 1:n
	for (int i = 0; i < n; i++) {
		
		// inInv(iIn:iIn+in(i) - 1) = i;
		VECTOR_T* inInv_ind = sequence(iIn, iIn + (int)VECTOR_ID(get)(in, i) - 1);
		if (inInv_ind != NULL) {
			ordinal_index_assign(inInv, inInv_ind, (FP_T)i);
			VECTOR_ID(free)(inInv_ind);
		}
		
		// outInv(iOut:iOut+out(i) - 1) = i;
		VECTOR_T* outInv_ind = sequence(iOut, iOut + (int)VECTOR_ID(get)(out, i) - 1);
		if (outInv_ind != NULL) {
			ordinal_index_assign(outInv, outInv_ind, (FP_T)i);
			VECTOR_ID(free)(outInv_ind);
		}
		
		// iIn = iIn+in(i);
		iIn += (int)VECTOR_ID(get)(in, i);
		
		// iOut = iOut+out(i);
		iOut += (int)VECTOR_ID(get)(out, i);
	}
	
	// cij = eye(n);
	MATRIX_T* cij = eye(n);
	
	// edges = [outInv(1:k)'; inInv(randperm(k))'];
	VECTOR_T* outInv_ind = sequence(0, k - 1);
	VECTOR_T* edges_row_0 = ordinal_index(outInv, outInv_ind);
	VECTOR_ID(free)(outInv);
	VECTOR_ID(free)(outInv_ind);
	gsl_permutation* inInv_ind = randperm(k);
	VECTOR_T* edges_row_1 = permute(inInv_ind, inInv);
	gsl_permutation_free(inInv_ind);
	VECTOR_ID(free)(inInv);
	MATRIX_T* edges = concatenate_columns(edges_row_0, edges_row_1);
	VECTOR_ID(free)(edges_row_0);
	VECTOR_ID(free)(edges_row_1);
	
	bool flag = true;
	
	// for i = 1:k
	for (int i = 0; i < k && flag; i++) {
		
		// if cij(edges(1,i),edges(2,i)),
		if (fp_nonzero(MATRIX_ID(get)(cij, (int)MATRIX_ID(get)(edges, 0, i), (int)MATRIX_ID(get)(edges, 1, i)))) {
			
			// warningCounter = 1;
			int warningCounter = 1;
			
			// while (1)
			while (true) {
				
				// switchTo = ceil(k*rand);
				int switchTo = (int)std::ceil((k - 1) * gsl_rng_uniform(rng));
				
				// if ~(cij(edges(1,i),edges(2,switchTo)) || cij(edges(1,switchTo),edges(2,i))),
				if (!(fp_nonzero(MATRIX_ID(get)(cij, (int)MATRIX_ID(get)(edges, 0, i), (int)MATRIX_ID(get)(edges, 1, switchTo))) ||
					  fp_nonzero(MATRIX_ID(get)(cij, (int)MATRIX_ID(get)(edges, 0, switchTo), (int)MATRIX_ID(get)(edges, 1, i))))) {
					
					// cij(edges(1,i),edges(2,switchTo)) = 1;
					MATRIX_ID(set)(cij, (int)MATRIX_ID(get)(edges, 0, i), (int)MATRIX_ID(get)(edges, 1, switchTo), 1.0);
					
					// if switchTo < i,
					if (switchTo < i) {
						
						// cij(edges(1,switchTo),edges(2,switchTo)) = 0;
						MATRIX_ID(set)(cij, (int)MATRIX_ID(get)(edges, 0, switchTo), (int)MATRIX_ID(get)(edges, 1, switchTo), 0.0);
						
						// cij(edges(1,switchTo),edges(2,i)) = 1;
						MATRIX_ID(set)(cij, (int)MATRIX_ID(get)(edges, 0, switchTo), (int)MATRIX_ID(get)(edges, 1, i), 1.0);
					}
					
					// temp = edges(2,i);
					FP_T temp = MATRIX_ID(get)(edges, 1, i);
					
					// edges(2,i) = edges(2,switchTo);
					MATRIX_ID(set)(edges, 1, i, MATRIX_ID(get)(edges, 1, switchTo));
					
					// edges(2,switchTo) = temp;
					MATRIX_ID(set)(edges, 1, switchTo, temp);
					
					// break
					break;
				}
				
				// warningCounter = warningCounter+1;
				warningCounter++;
				
				// if warningCounter == 2*k^2
				if (warningCounter == 2 * k * k) {
					
					// flag = 0;
					flag = false;
					
					// return;
					break;
				}
			}
		} else {
			
			// cij(edges(1,i),edges(2,i)) = 1;
			MATRIX_ID(set)(cij, (int)MATRIX_ID(get)(edges, 0, i), (int)MATRIX_ID(get)(edges, 1, i), 1.0);
		}
	}

	MATRIX_ID(free)(edges);
	
	// flag = 1;
	if (!flag) {
		MATRIX_ID(free)(cij);
		return NULL;
	}
	
	// cij = cij - eye(n);
	MATRIX_T* eye_n = eye(n);
	MATRIX_ID(sub)(cij, eye_n);
	MATRIX_ID(free)(eye_n);
	
	return cij;
}

/*
 * Generates a random directed binary graph with the same in-degree and out-
 * degree sequences of the given graph.  Since the degree sequences are
 * necessarily valid, this function should not return NULL unless the given
 * matrix contains nonzero entries on the main diagonal.
 */
MATRIX_T* BCT_NAMESPACE::makerandCIJdegreesfixed(const MATRIX_T* m) {
	if (safe_mode) check_status(m, SQUARE | NO_LOOPS, "makerandCIJdegreesfixed");
	MATRIX_T* ret;
	do {
		VECTOR_T* id;
		VECTOR_T* od;
		VECTOR_T* deg = degrees_dir(m, &id, &od);
		VECTOR_ID(free)(deg);
		ret = makerandCIJdegreesfixed(id, od);
		VECTOR_ID(free)(id);
		VECTOR_ID(free)(od);
	} while (ret == NULL && has_no_loops(m));
	return ret;
}
