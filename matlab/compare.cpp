#include <cmath>

#include "matlab.h"

FP_T MATLAB_NAMESPACE::epsilon = std::numeric_limits<FP_T>::epsilon();

/*
 * Compares two floating-point numbers.
 */
int MATLAB_NAMESPACE::fp_compare(FP_T x, FP_T y) {
	if (fp_zero(x) && fp_zero(y)) {
		return 0;
	} else {
		int exponent;
		FP_T max = (std::abs(x) > std::abs(y)) ? x : y;
		std::frexp(max, &exponent);
		FP_T delta = std::ldexp(epsilon, exponent);
		FP_T difference = x - y;
		if (difference > delta) {
			return 1;
		} else if (difference < -delta) {
			return -1;
		} else {
			return 0;
		}
	}
}

bool MATLAB_NAMESPACE::fp_zero(FP_T x) { return std::abs(x) < epsilon; }
bool MATLAB_NAMESPACE::fp_nonzero(FP_T x) { return std::abs(x) > epsilon; }
bool MATLAB_NAMESPACE::fp_equal(FP_T x, FP_T y) { return fp_compare(x, y) == 0; }
bool MATLAB_NAMESPACE::fp_not_equal(FP_T x, FP_T y) { return fp_compare(x, y) != 0; }
bool MATLAB_NAMESPACE::fp_less(FP_T x, FP_T y) { return fp_compare(x, y) == -1; }
bool MATLAB_NAMESPACE::fp_less_or_equal(FP_T x, FP_T y) { return fp_compare(x, y) <= 0; }
bool MATLAB_NAMESPACE::fp_greater(FP_T x, FP_T y) { return fp_compare(x, y) == 1; }
bool MATLAB_NAMESPACE::fp_greater_or_equal(FP_T x, FP_T y) { return fp_compare(x, y) >= 0; }

/*
 * Compares two vectors lexicographically, returning -1, 0, or 1 if the first
 * vector is less than, equal to, or greater than the second.
 */
int MATLAB_NAMESPACE::compare_vectors(const VECTOR_T* v1, const VECTOR_T* v2) {
	for (int i = 0; i < (int)v1->size; i++) {
		if (i >= (int)v2->size) {
			return 1;
		}
		int result = fp_compare(VECTOR_ID(get)(v1, i), VECTOR_ID(get)(v2, i));
		if (result != 0) {
			return result;
		}
	}
	if (v1->size < v2->size) {
		return -1;
	} else {
		return 0;
	}
}

/*
 * Returns whether the first vector comes before the second vector in a strict
 * weak ordering.
 */
bool MATLAB_NAMESPACE::vector_less(VECTOR_T* v1, VECTOR_T* v2) {
	return compare_vectors(v1, v2) == -1;
}

/*
 * Compares two matrices lexicographically, returning -1, 0, or 1 if the first
 * matrix is less than, equal to, or greater than the second.
 */
int MATLAB_NAMESPACE::compare_matrices(const MATRIX_T* m1, const MATRIX_T* m2) {
	int size1 = (int)m1->size1 * (int)m1->size2;
	int size2 = (int)m2->size1 * (int)m2->size2;
	for (int i = 0; i < size1; i++) {
		if (i >= size2) {
			return 1;
		}
		int result = fp_compare(ordinal_index(m1, i), ordinal_index(m2, i));
		if (result != 0) {
			return result;
		}
	}
	if (size1 < size2) {
		return -1;
	} else {
		return 0;
	}
}

/*
 * Returns whether the first matrix comes before the second matrix in a strict
 * weak ordering.
 */
bool MATLAB_NAMESPACE::matrix_less(MATRIX_T* m1, MATRIX_T* m2) {
	return compare_matrices(m1, m2) == -1;
}

/*
 * Emulates (v op x), where op is a binary comparison operator.
 */
VECTOR_T* MATLAB_NAMESPACE::compare_elements(const VECTOR_T* v, comparator compare, FP_T x) {
	VECTOR_T* cmp_v = VECTOR_ID(alloc)(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		FP_T value = VECTOR_ID(get)(v, i);
		VECTOR_ID(set)(cmp_v, i, (FP_T)compare(value, x));
	}
	return cmp_v;
}

/*
 * Emulates (v1 op v2), where op is a binary comparison operator.
 */
VECTOR_T* MATLAB_NAMESPACE::compare_elements(const VECTOR_T* v1, comparator compare, const VECTOR_T* v2) {
	if (v1->size != v2->size) {
		return NULL;
	}
	VECTOR_T* cmp_v = VECTOR_ID(alloc)(v1->size);
	for (int i = 0; i < (int)v1->size; i++) {
		FP_T value1 = VECTOR_ID(get)(v1, i);
		FP_T value2 = VECTOR_ID(get)(v2, i);
		VECTOR_ID(set)(cmp_v, i, (FP_T)compare(value1, value2));
	}
	return cmp_v;
}

/*
 * Emulates (m op x), where op is a binary comparison operator.
 */
MATRIX_T* MATLAB_NAMESPACE::compare_elements(const MATRIX_T* m, comparator compare, FP_T x) {
	MATRIX_T* cmp_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			FP_T value = MATRIX_ID(get)(m, i, j);
			MATRIX_ID(set)(cmp_m, i, j, (FP_T)compare(value, x));
		}
	}
	return cmp_m;
}

/*
 * Emulates (m1 op m2), where op is a binary comparison operator.
 */
MATRIX_T* MATLAB_NAMESPACE::compare_elements(const MATRIX_T* m1, comparator compare, const MATRIX_T* m2) {
	if (m1->size1 != m2->size1 || m1->size2 != m2->size2) {
		return NULL;
	}
	MATRIX_T* cmp_m = MATRIX_ID(alloc)(m1->size1, m1->size2);
	for (int i = 0; i < (int)m1->size1; i++) {
		for (int j = 0; j < (int)m1->size2; j++) {
			FP_T value1 = MATRIX_ID(get)(m1, i, j);
			FP_T value2 = MATRIX_ID(get)(m2, i, j);
			MATRIX_ID(set)(cmp_m, i, j, (FP_T)compare(value1, value2));
		}
	}
	return cmp_m;
}
