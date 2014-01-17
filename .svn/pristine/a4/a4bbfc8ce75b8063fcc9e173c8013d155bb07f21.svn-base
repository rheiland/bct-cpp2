#include <cmath>

#include "bct.h"

/*
 * Returns a copy of the given matrix with each nonzero element inverted.
 */
MATRIX_T* BCT_NAMESPACE::invert_elements(const MATRIX_T* m) {
	MATRIX_T* inv_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			FP_T value = MATRIX_ID(get)(m, i, j);
			if (fp_nonzero(value)) {
				MATRIX_ID(set)(inv_m, i, j, 1.0 / value);
			} else {
				MATRIX_ID(set)(inv_m, i, j, 0.0);
			}
		}
	}
	return inv_m;
}

/*
 * Returns a copy of the given matrix with no loops.
 */
MATRIX_T* BCT_NAMESPACE::remove_loops(const MATRIX_T* m) {
	MATRIX_T* nl_m = copy(m);
	VECTOR_ID(view) diag_nl_m = MATRIX_ID(diagonal)(nl_m);
	VECTOR_ID(set_zero)(&diag_nl_m.vector);
	return nl_m;
}

/*
 * Returns a binary copy of the given matrix.
 */
MATRIX_T* BCT_NAMESPACE::to_binary(const MATRIX_T* m) {
	return compare_elements(m, fp_not_equal, 0.0);
}

/*
 * Returns a positive copy of the given matrix.
 */
MATRIX_T* BCT_NAMESPACE::to_positive(const MATRIX_T* m) {
	MATRIX_T* pos_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			MATRIX_ID(set)(pos_m, i, j, std::abs(MATRIX_ID(get)(m, i, j)));
		}
	}
	return pos_m;
}

/*
 * Returns an undirected copy of the given binary matrix.  For every pair of
 * nodes, if either m(i, j) or m(j, i) is nonzero, then both m(i, j) and m(j, i)
 * are set to one.  Otherwise, both are set to zero.
 */
MATRIX_T* BCT_NAMESPACE::to_undirected_bin(const MATRIX_T* m) {
	MATRIX_T* und_m = MATRIX_ID(calloc)(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = i; j < (int)m->size2; j++) {
			FP_T value_ij = MATRIX_ID(get)(m, i, j);
			FP_T value_ji = MATRIX_ID(get)(m, j, i);
			if (fp_nonzero(value_ij) || fp_nonzero(value_ji)) {
				MATRIX_ID(set)(und_m, i, j, 1.0);
				MATRIX_ID(set)(und_m, j, i, 1.0);
			}
		}
	}
	return und_m;
}

/*
 * Returns an undirected copy of the given weighted matrix.  For every pair of
 * nodes, m(i, j) and m(j, i) are both set to the average of their two values.
 */
MATRIX_T* BCT_NAMESPACE::to_undirected_wei(const MATRIX_T* m) {
	MATRIX_T* und_m = MATRIX_ID(calloc)(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = i; j < (int)m->size2; j++) {
			FP_T value_ij = MATRIX_ID(get)(m, i, j);
			FP_T value_ji = MATRIX_ID(get)(m, j, i);
			if (fp_nonzero(value_ij) || fp_nonzero(value_ji)) {
				FP_T average = (value_ij + value_ji) / 2.0;
				MATRIX_ID(set)(und_m, i, j, average);
				MATRIX_ID(set)(und_m, j, i, average);
			}
		}
	}
	return und_m;
}
