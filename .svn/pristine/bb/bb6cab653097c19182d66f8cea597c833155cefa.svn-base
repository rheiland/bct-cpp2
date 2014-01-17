#include <ctime>

#include "matlab.h"

/*
 * Returns a random number generator that is guaranteed to be seeded only once
 * during program execution.  This generator should not be freed by the caller.
 */
gsl_rng* MATLAB_NAMESPACE::get_rng() {
	static gsl_rng* rng = NULL;
	if (rng == NULL) {
		gsl_rng_default_seed = std::time(NULL);
		rng = gsl_rng_alloc(gsl_rng_default);
	}
	return rng;
}

/*
 * Seeds the given random number generator.
 */
void MATLAB_NAMESPACE::seed_rng(const gsl_rng* rng, unsigned long seed) {
	gsl_rng_set(rng, seed);
}

/*
 * Permutes the elements of a vector.
 */
VECTOR_T* MATLAB_NAMESPACE::permute(const gsl_permutation* p, const VECTOR_T* v) {
	if (p->size != v->size) return NULL;
	VECTOR_T* permuted_v = VECTOR_ID(alloc)(v->size);
	for (int i = 0; i < (int)p->size; i++) {
		int index = gsl_permutation_get(p, i);
		FP_T value = VECTOR_ID(get)(v, index);
		VECTOR_ID(set)(permuted_v, i, value);
	}
	return permuted_v;
}

/*
 * Permutes the columns of a matrix.
 */
MATRIX_T* MATLAB_NAMESPACE::permute_columns(const gsl_permutation* p, const MATRIX_T* m) {
	if (p->size != m->size2) return NULL;
	MATRIX_T* permuted_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0; i < (int)p->size; i++) {
		int i_col = gsl_permutation_get(p, i);
		VECTOR_ID(const_view) m_col_i_col = MATRIX_ID(const_column)(m, i_col);
		MATRIX_ID(set_col)(permuted_m, i, &m_col_i_col.vector);
	}
	return permuted_m;
}

/*
 * Permutes the rows of a matrix.
 */
MATRIX_T* MATLAB_NAMESPACE::permute_rows(const gsl_permutation* p, const MATRIX_T* m) {
	if (p->size != m->size1) return NULL;
	MATRIX_T* permuted_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0; i < (int)p->size; i++) {
		int i_row = gsl_permutation_get(p, i);
		VECTOR_ID(const_view) m_row_i_row = MATRIX_ID(const_row)(m, i_row);
		MATRIX_ID(set_row)(permuted_m, i, &m_row_i_row.vector);
	}
	return permuted_m;
}
