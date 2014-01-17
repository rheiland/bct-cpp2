#include <cmath>
#include <gsl/gsl_blas.h>

#include "matlab.h"

/*
 * Emulates ([v x]) for a row vector or ([v ; x]) for a column vector.
 */
VECTOR_T* MATLAB_NAMESPACE::concatenate(const VECTOR_T* v, FP_T x) {
	if (v == NULL) {
		VECTOR_T* cat_v = VECTOR_ID(alloc)(1);
		VECTOR_ID(set)(cat_v, 0, x);
		return cat_v;
	}
	VECTOR_T* cat_v = VECTOR_ID(alloc)(v->size + 1);
	VECTOR_ID(view) cat_subv = VECTOR_ID(subvector)(cat_v, 0, v->size);
	VECTOR_ID(memcpy)(&cat_subv.vector, v);
	VECTOR_ID(set)(cat_v, v->size, x);
	return cat_v;
}

/*
 * Emulates ([x v]) for a row vector or ([x ; v]) for a column vector.
 */
VECTOR_T* MATLAB_NAMESPACE::concatenate(FP_T x, const VECTOR_T* v) {
	if (v == NULL) {
		VECTOR_T* cat_v = VECTOR_ID(alloc)(1);
		VECTOR_ID(set)(cat_v, 0, x);
		return cat_v;
	}
	VECTOR_T* cat_v = VECTOR_ID(alloc)(v->size + 1);
	VECTOR_ID(view) cat_subv = VECTOR_ID(subvector)(cat_v, 1, v->size);
	VECTOR_ID(memcpy)(&cat_subv.vector, v);
	VECTOR_ID(set)(cat_v, 0, x);
	return cat_v;
}
 
/*
 * Emulates ([v1 v2]) for row vectors or ([v1 ; v2]) for column vectors.
 */
VECTOR_T* MATLAB_NAMESPACE::concatenate(const VECTOR_T* v1, const VECTOR_T* v2) {
	if (v1 == NULL && v2 == NULL) {
		return NULL;
	} else if (v1 == NULL) {
		return copy(v2);
	} else if (v2 == NULL) {
		return copy(v1);
	}
	VECTOR_T* cat_v = VECTOR_ID(alloc)(v1->size + v2->size);
	VECTOR_ID(view) cat_subv1 = VECTOR_ID(subvector)(cat_v, 0, v1->size);
	VECTOR_ID(view) cat_subv2 = VECTOR_ID(subvector)(cat_v, v1->size, v2->size);
	VECTOR_ID(memcpy)(&cat_subv1.vector, v1);
	VECTOR_ID(memcpy)(&cat_subv2.vector, v2);
	return cat_v;
}

/*
 * Emulates ([v1 ; v2]) for row vectors.
 */
MATRIX_T* MATLAB_NAMESPACE::concatenate_columns(const VECTOR_T* v1, const VECTOR_T* v2) {
	if (v1 == NULL && v2 == NULL) {
		return NULL;
	} else if (v1 == NULL) {
		return to_row_matrix(v2);
	} else if (v2 == NULL) {
		return to_row_matrix(v1);
	} else if (v1->size != v2->size) {
		return NULL;
	}
	MATRIX_T* cat_m = MATRIX_ID(alloc)(2, v1->size);
	MATRIX_ID(set_row)(cat_m, 0, v1);
	MATRIX_ID(set_row)(cat_m, 1, v2);
	return cat_m;
}

/*
 * Emulates ([m ; v]) for a row vector.
 */
MATRIX_T* MATLAB_NAMESPACE::concatenate_columns(const MATRIX_T* m, const VECTOR_T* v) {
	if (m == NULL && v == NULL) {
		return NULL;
	} else if (m == NULL) {
		return to_row_matrix(v);
	} else if (v == NULL) {
		return copy(m);
	} else if (m->size2 != v->size) {
		return NULL;
	}
	MATRIX_T* cat_m = MATRIX_ID(alloc)(m->size1 + 1, m->size2);
	MATRIX_ID(view) cat_subm = MATRIX_ID(submatrix)(cat_m, 0, 0, m->size1, m->size2);
	MATRIX_ID(memcpy)(&cat_subm.matrix, m);
	MATRIX_ID(set_row)(cat_m, m->size1, v);
	return cat_m;
}

/*
 * Emulates ([v ; m]) for a row vector.
 */
MATRIX_T* MATLAB_NAMESPACE::concatenate_columns(const VECTOR_T* v, const MATRIX_T* m) {
	if (m == NULL && v == NULL) {
		return NULL;
	} else if (m == NULL) {
		return to_row_matrix(v);
	} else if (v == NULL) {
		return copy(m);
	} else if (m->size2 != v->size) {
		return NULL;
	}
	MATRIX_T* cat_m = MATRIX_ID(alloc)(m->size1 + 1, m->size2);
	MATRIX_ID(set_row)(cat_m, 0, v);
	MATRIX_ID(view) cat_subm = MATRIX_ID(submatrix)(cat_m, 1, 0, m->size1, m->size2);
	MATRIX_ID(memcpy)(&cat_subm.matrix, m);
	return cat_m;
}

/*
 * Emulates ([m1 ; m2]).
 */
MATRIX_T* MATLAB_NAMESPACE::concatenate_columns(const MATRIX_T* m1, const MATRIX_T* m2) {
	if (m1 == NULL && m2 == NULL) {
		return NULL;
	} else if (m1 == NULL) {
		return copy(m2);
	} else if (m2 == NULL) {
		return copy(m1);
	} else if (m1->size2 != m2->size2) {
		return NULL;
	}
	MATRIX_T* cat_m = MATRIX_ID(alloc)(m1->size1 + m2->size1, m1->size2);
	MATRIX_ID(view) cat_subm1 = MATRIX_ID(submatrix)(cat_m, 0, 0, m1->size1, m1->size2);
	MATRIX_ID(view) cat_subm2 = MATRIX_ID(submatrix)(cat_m, m1->size1, 0, m2->size1, m2->size2);
	MATRIX_ID(memcpy)(&cat_subm1.matrix, m1);
	MATRIX_ID(memcpy)(&cat_subm2.matrix, m2);
	return cat_m;
}

/*
 * Emulates ([v1 v2]) for column vectors.
 */
MATRIX_T* MATLAB_NAMESPACE::concatenate_rows(const VECTOR_T* v1, const VECTOR_T* v2) {
	if (v1 == NULL && v2 == NULL) {
		return NULL;
	} else if (v1 == NULL) {
		return to_column_matrix(v2);
	} else if (v2 == NULL) {
		return to_column_matrix(v1);
	} else if (v1->size != v2->size) {
		return NULL;
	}
	MATRIX_T* cat_m = MATRIX_ID(alloc)(v1->size, 2);
	MATRIX_ID(set_col)(cat_m, 0, v1);
	MATRIX_ID(set_col)(cat_m, 1, v2);
	return cat_m;
}

/*
 * Emulates ([m v]) for a column vector.
 */
MATRIX_T* MATLAB_NAMESPACE::concatenate_rows(const MATRIX_T* m, const VECTOR_T* v) {
	if (m == NULL && v == NULL) {
		return NULL;
	} else if (m == NULL) {
		return to_column_matrix(v);
	} else if (v == NULL) {
		return copy(m);
	} else if (m->size1 != v->size) {
		return NULL;
	}
	MATRIX_T* cat_m = MATRIX_ID(alloc)(m->size1, m->size2 + 1);
	MATRIX_ID(view) cat_subm = MATRIX_ID(submatrix)(cat_m, 0, 0, m->size1, m->size2);
	MATRIX_ID(memcpy)(&cat_subm.matrix, m);
	MATRIX_ID(set_col)(cat_m, m->size2, v);
	return cat_m;
}

/*
 * Emulates ([v m]) for a column vector.
 */
MATRIX_T* MATLAB_NAMESPACE::concatenate_rows(const VECTOR_T* v, const MATRIX_T* m) {
	if (m == NULL && v == NULL) {
		return NULL;
	} else if (m == NULL) {
		return to_column_matrix(v);
	} else if (v == NULL) {
		return copy(m);
	} else if (m->size1 != v->size) {
		return NULL;
	}
	MATRIX_T* cat_m = MATRIX_ID(alloc)(m->size1, m->size2 + 1);
	MATRIX_ID(set_col)(cat_m, 0, v);
	MATRIX_ID(view) cat_subm = MATRIX_ID(submatrix)(cat_m, 0, 1, m->size1, m->size2);
	MATRIX_ID(memcpy)(&cat_subm.matrix, m);
	return cat_m;
}

/*
 * Emulates ([m1 m2]).
 */
MATRIX_T* MATLAB_NAMESPACE::concatenate_rows(const MATRIX_T* m1, const MATRIX_T* m2) {
	if (m1 == NULL && m2 == NULL) {
		return NULL;
	} else if (m1 == NULL) {
		return copy(m2);
	} else if (m2 == NULL) {
		return copy(m1);
	} else if (m1->size1 != m2->size1) {
		return NULL;
	}
	MATRIX_T* cat_m = MATRIX_ID(alloc)(m1->size1, m1->size2 + m2->size2);
	MATRIX_ID(view) cat_subm1 = MATRIX_ID(submatrix)(cat_m, 0, 0, m1->size1, m1->size2);
	MATRIX_ID(view) cat_subm2 = MATRIX_ID(submatrix)(cat_m, 0, m1->size2, m2->size1, m2->size2);
	MATRIX_ID(memcpy)(&cat_subm1.matrix, m1);
	MATRIX_ID(memcpy)(&cat_subm2.matrix, m2);
	return cat_m;
}

/*
 * Emulates copy assignment.
 */
VECTOR_T* MATLAB_NAMESPACE::copy(const VECTOR_T* v) {
	VECTOR_T* copy_v = VECTOR_ID(alloc)(v->size);
	VECTOR_ID(memcpy)(copy_v, v);
	return copy_v;
}

/*
 * Emulates copy assignment.
 */
MATRIX_T* MATLAB_NAMESPACE::copy(const MATRIX_T* m) {
	MATRIX_T* copy_m = MATRIX_ID(alloc)(m->size1, m->size2);
	MATRIX_ID(memcpy)(copy_m, m);
	return copy_m;
}

/*
 * Emulates (m1 \ m2) = (inv(m1) * m2).
 */
MATRIX_T* MATLAB_NAMESPACE::div_left(const MATRIX_T* m1, const MATRIX_T* m2) {
	if (m1->size1 != m1->size2 || m2->size1 != m2->size2 || m1->size1 != m2->size1) {
		return NULL;
	}
	MATRIX_T* inv_m1 = inv(m1);
	MATRIX_T* div_m = mul(inv_m1, m2);
	MATRIX_ID(free)(inv_m1);
	return div_m;
}

/*
 * Emulates (m1 / m2) = ((inv(m2') * m1')').
 */
MATRIX_T* MATLAB_NAMESPACE::div_right(const MATRIX_T* m1, const MATRIX_T* m2) {
	if (m1->size1 != m1->size2 || m2->size1 != m2->size2 || m1->size1 != m2->size1) {
		return NULL;
	}
	MATRIX_T* m2_transpose = MATRIX_ID(alloc)(m2->size2, m2->size1);
	MATRIX_ID(transpose_memcpy)(m2_transpose, m2);
	MATRIX_T* inv_m2_transpose = inv(m2_transpose);
	MATRIX_ID(free)(m2_transpose);
	MATRIX_T* m1_transpose = MATRIX_ID(alloc)(m1->size2, m1->size1);
	MATRIX_ID(transpose_memcpy)(m1_transpose, m1);
	MATRIX_T* div_m = mul(inv_m2_transpose, m1_transpose);
	MATRIX_ID(free)(inv_m2_transpose);
	MATRIX_ID(free)(m1_transpose);
	MATRIX_ID(transpose)(div_m);
	return div_m;
}

/*
 * Emulates (v1 & v2).
 */
VECTOR_T* MATLAB_NAMESPACE::logical_and(const VECTOR_T* v1, const VECTOR_T* v2) {
	if (v1->size != v2->size) {
		return NULL;
	}
	VECTOR_T* and_v = VECTOR_ID(alloc)(v1->size);
	for (int i = 0; i < (int)v1->size; i++) {
		bool nz1 = fp_nonzero(VECTOR_ID(get)(v1, i));
		bool nz2 = fp_nonzero(VECTOR_ID(get)(v2, i));
		VECTOR_ID(set)(and_v, i, (FP_T)(nz1 && nz2));
	}
	return and_v;
}

/*
 * Emulates (m1 & m2).
 */
MATRIX_T* MATLAB_NAMESPACE::logical_and(const MATRIX_T* m1, const MATRIX_T* m2) {
	if (m1->size1 != m2->size1 || m1->size2 != m2->size2) {
		return NULL;
	}
	MATRIX_T* and_m = MATRIX_ID(alloc)(m1->size1, m1->size2);
	for (int i = 0; i < (int)m1->size1; i++) {
		for (int j = 0; j < (int)m1->size2; j++) {
			bool nz1 = fp_nonzero(MATRIX_ID(get)(m1, i, j));
			bool nz2 = fp_nonzero(MATRIX_ID(get)(m2, i, j));
			MATRIX_ID(set)(and_m, i, j, (FP_T)(nz1 && nz2));
		}
	}
	return and_m;
}

/*
 * Emulates (~v).
 */
VECTOR_T* MATLAB_NAMESPACE::logical_not(const VECTOR_T* v) {
	VECTOR_T* not_v = VECTOR_ID(alloc)(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		bool z = fp_zero(VECTOR_ID(get)(v, i));
		VECTOR_ID(set)(not_v, i, (FP_T)z);
	}
	return not_v;
}

/*
 * Emulates (~m)
 */
MATRIX_T* MATLAB_NAMESPACE::logical_not(const MATRIX_T* m) {
	MATRIX_T* not_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			bool z = fp_zero(MATRIX_ID(get)(m, i, j));
			MATRIX_ID(set)(not_m, i, j, (FP_T)z);
		}
	}
	return not_m;
}

/*
 * Emulates (v1 | v2).
 */
VECTOR_T* MATLAB_NAMESPACE::logical_or(const VECTOR_T* v1, const VECTOR_T* v2) {
	if (v1->size != v2->size) {
		return NULL;
	}
	VECTOR_T* or_v = VECTOR_ID(alloc)(v1->size);
	for (int i = 0; i < (int)v1->size; i++) {
		bool nz1 = fp_nonzero(VECTOR_ID(get)(v1, i));
		bool nz2 = fp_nonzero(VECTOR_ID(get)(v2, i));
		VECTOR_ID(set)(or_v, i, (FP_T)(nz1 || nz2));
	}
	return or_v;
}

/*
 * Emulates (m1 | m2).
 */
MATRIX_T* MATLAB_NAMESPACE::logical_or(const MATRIX_T* m1, const MATRIX_T* m2) {
	if (m1->size1 != m2->size1 || m1->size2 != m2->size2) {
		return NULL;
	}
	MATRIX_T* or_m = MATRIX_ID(alloc)(m1->size1, m1->size2);
	for (int i = 0; i < (int)m1->size1; i++) {
		for (int j = 0; j < (int)m1->size2; j++) {
			bool nz1 = fp_nonzero(MATRIX_ID(get)(m1, i, j));
			bool nz2 = fp_nonzero(MATRIX_ID(get)(m2, i, j));
			MATRIX_ID(set)(or_m, i, j, (FP_T)(nz1 || nz2));
		}
	}
	return or_m;
}

/*
 * Emulates (m1 * m2).
 */
MATRIX_T* MATLAB_NAMESPACE::mul(const MATRIX_T* m1, const MATRIX_T* m2) {
	if (m1->size2 != m2->size1) {
		return NULL;
	}
	MATRIX_T* mul_m = MATRIX_ID(alloc)(m1->size1, m2->size2);
	for (int i = 0; i < (int)m1->size1; i++) {
		VECTOR_ID(const_view) row = MATRIX_ID(const_row)(m1, i);
		for (int j = 0; j < (int)m2->size2; j++) {
			VECTOR_ID(const_view) col = MATRIX_ID(const_column)(m2, j);
			FP_T value = 0.0;
			for (int k = 0; k < (int)m1->size2; k++) {
				value += VECTOR_ID(get)(&row.vector, k) * VECTOR_ID(get)(&col.vector, k);
			}
			MATRIX_ID(set)(mul_m, i, j, value);
		}
	}
	return mul_m;
}

/*
 * Emulates (m ^ power).
 */
MATRIX_T* MATLAB_NAMESPACE::pow(const MATRIX_T* m, int power) {
	if (m->size1 != m->size2 || power < 1) {
		return NULL;
	}
	MATRIX_T* pow_m = copy(m);
	for (int i = 2; i <= power; i++) {
		MATRIX_T* temp_m = mul(pow_m, m);
		MATRIX_ID(free)(pow_m);
		pow_m = temp_m;
	}
	return pow_m;
}

/*
 * Emulates (v .^ power).
 */
VECTOR_T* MATLAB_NAMESPACE::pow_elements(const VECTOR_T* v, FP_T power) {
	VECTOR_T* pow_v = VECTOR_ID(alloc)(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		FP_T value = std::pow(VECTOR_ID(get)(v, i), power);
		VECTOR_ID(set)(pow_v, i, value);
	}
	return pow_v;
}

/*
 * Emulates (v .^ powers).
 */
VECTOR_T* MATLAB_NAMESPACE::pow_elements(const VECTOR_T* v, const VECTOR_T* powers) {
	if (v->size != powers->size) {
		return NULL;
	}
	VECTOR_T* pow_v = VECTOR_ID(alloc)(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		FP_T value = std::pow(VECTOR_ID(get)(v, i), VECTOR_ID(get)(powers, i));
		VECTOR_ID(set)(pow_v, i, value);
	}
	return pow_v;
}

/*
 * Emulates (m .^ power).
 */
MATRIX_T* MATLAB_NAMESPACE::pow_elements(const MATRIX_T* m, FP_T power) {
	MATRIX_T* pow_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			FP_T value = std::pow(MATRIX_ID(get)(m, i, j), power);
			MATRIX_ID(set)(pow_m, i, j, value);
		}
	}
	return pow_m;
}

/*
 * Emulates (m .^ powers).
 */
MATRIX_T* MATLAB_NAMESPACE::pow_elements(const MATRIX_T* m, const MATRIX_T* powers) {
	if (m->size1 != powers->size1 || m->size2 != powers->size2) {
		return NULL;
	}
	MATRIX_T* pow_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			FP_T value = std::pow(MATRIX_ID(get)(m, i, j), MATRIX_ID(get)(powers, i, j));
			MATRIX_ID(set)(pow_m, i, j, value);
		}
	}
	return pow_m;
}

/* 
 * Emulates (start:end).
 */
VECTOR_T* MATLAB_NAMESPACE::sequence(int start, int end) {
	return sequence(start, 1, end);
}

/*
 * Emulates (start:step:end).
 */
VECTOR_T* MATLAB_NAMESPACE::sequence(int start, int step, int end) {
	int n_seq = (end - start) / step + 1;
	if (n_seq <= 0) {
		return NULL;
	}
	VECTOR_T* seq_v = VECTOR_ID(alloc)(n_seq);
	for (int i = 0, value = start; i < n_seq; i++, value += step) {
		VECTOR_ID(set)(seq_v, i, value);
	}
	return seq_v;
}
