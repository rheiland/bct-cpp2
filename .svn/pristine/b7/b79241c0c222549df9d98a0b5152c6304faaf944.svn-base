#include "matlab.h"

/*
 * Converts a vector to an array.
 */
void MATLAB_NAMESPACE::to_array(const VECTOR_T* v, FP_T* array) {
	for (int i = 0; i < (int)v->size; i++) {
		array[i] = VECTOR_ID(get)(v, i);
	}
}

/*
 * Converts a vector to a boolean: true if all elements are nonzero, false
 * otherwise.
 */
bool MATLAB_NAMESPACE::to_bool(const VECTOR_T* v) {
	return all(v) == 1;
}

/*
 * Converts a matrix to a boolean: true if all elements are nonzero, false
 * otherwise.
 */
bool MATLAB_NAMESPACE::to_bool(const MATRIX_T* m) {
	VECTOR_T* all_v = all(m);
	bool ret = all(all_v);
	VECTOR_ID(free)(all_v);
	return ret;
}

/*
 * Converts a matrix to a vector.  The vector is constructed by consecutively
 * appending columns.
 */
VECTOR_T* MATLAB_NAMESPACE::to_vector(const MATRIX_T* m) {
	VECTOR_T* v = VECTOR_ID(alloc)(m->size1 * m->size2);
	for (int j = 0; j < (int)m->size2; j++) {
		for (int i = 0; i < (int)m->size1; i++) {
			FP_T value = MATRIX_ID(get)(m, i, j);
			VECTOR_ID(set)(v, j * m->size1 + i, value);
		}
	}
	return v;
}

/*
 * Converts a double-precision vector to the currently selected precision.
 */
VECTOR_T* MATLAB_NAMESPACE::to_vector(const gsl_vector* v_d) {
	VECTOR_T* v = VECTOR_ID(alloc)(v_d->size);
	for (int i = 0; i < (int)v_d->size; i++) {
		FP_T value = (FP_T)gsl_vector_get(v_d, i);
		VECTOR_ID(set)(v, i, value);
	}
	return v;
}

/*
 * Converts a vector to double precision.
 */
gsl_vector* MATLAB_NAMESPACE::to_vector_double(const VECTOR_T* v) {
	gsl_vector* v_d = gsl_vector_alloc(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		double value = (double)VECTOR_ID(get)(v, i);
		gsl_vector_set(v_d, i, value);
	}
	return v_d;
}

/*
 * Converts a vector to a single-column matrix.
 */
MATRIX_T* MATLAB_NAMESPACE::to_column_matrix(const VECTOR_T* v) {
	MATRIX_T* m = MATRIX_ID(alloc)(v->size, 1);
	for (int i = 0; i < (int)v->size; i++) {
		MATRIX_ID(set)(m, i, 0, VECTOR_ID(get)(v, i));
	}
	return m;
}

/*
 * Converts a vector to a single-row matrix.
 */
MATRIX_T* MATLAB_NAMESPACE::to_row_matrix(const VECTOR_T* v) {
	MATRIX_T* m = MATRIX_ID(alloc)(1, v->size);
	for (int i = 0; i < (int)v->size; i++) {
		MATRIX_ID(set)(m, 0, i, VECTOR_ID(get)(v, i));
	}
	return m;
}

/*
 * Converts a double-precision matrix to the currently selected precision.
 */
MATRIX_T* MATLAB_NAMESPACE::to_matrix(const gsl_matrix* m_d) {
	MATRIX_T* m = MATRIX_ID(alloc)(m_d->size1, m_d->size2);
	for (int i = 0; i < (int)m_d->size1; i++) {
		for (int j = 0; j < (int)m_d->size2; j++) {
			FP_T value = (FP_T)gsl_matrix_get(m_d, i, j);
			MATRIX_ID(set)(m, i, j, value);
		}
	}
	return m;
}

/*
 * Converts a matrix to double precision.
 */
gsl_matrix* MATLAB_NAMESPACE::to_matrix_double(const MATRIX_T* m) {
	gsl_matrix* m_d = gsl_matrix_alloc(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			double value = (double)MATRIX_ID(get)(m, i, j);
			gsl_matrix_set(m_d, i, j, value);
		}
	}
	return m_d;
}

/*
 * Converts a permutation to a vector.
 */
VECTOR_T* MATLAB_NAMESPACE::to_vector(const gsl_permutation* p) {
	VECTOR_T* v = VECTOR_ID(alloc)(p->size);
	for (int i = 0; i < (int)p->size; i++) {
		VECTOR_ID(set)(v, i, (FP_T)gsl_permutation_get(p, i));
	}
	return v;
}

/*
 * Converts a vector to a permutation.
 */
gsl_permutation* MATLAB_NAMESPACE::to_permutation(const VECTOR_T* v) {
	gsl_permutation* p = gsl_permutation_alloc(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		p->data[i] = (int)VECTOR_ID(get)(v, i);
	}
	if (gsl_permutation_valid(p) == 1) {
		gsl_permutation_free(p);
		return NULL;
	} else {
		return p;
	}
}
