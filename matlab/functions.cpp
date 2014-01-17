#include <cmath>
#include <cstddef>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "matlab.h"
#include "sort.h"

/*
 * See MATLAB documentation for descriptions of these functions.  We will
 * document instances where our version differs from the MATLAB version.
 */

VECTOR_T* MATLAB_NAMESPACE::abs(const VECTOR_T* v) {
	VECTOR_T* abs_v = VECTOR_ID(alloc)(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		VECTOR_ID(set)(abs_v, i, std::abs(VECTOR_ID(get)(v, i)));
	}
	return abs_v;
}

MATRIX_T* MATLAB_NAMESPACE::abs(const MATRIX_T* m) {
	MATRIX_T* abs_m = MATRIX_ID(alloc)(m->size1, m->size2);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			MATRIX_ID(set)(abs_m, i, j, std::abs(MATRIX_ID(get)(m, i, j)));
		}
	}
	return abs_m;
}

int MATLAB_NAMESPACE::all(const VECTOR_T* v) {
	for (int i = 0; i < (int)v->size; i++) {
		if (fp_zero(VECTOR_ID(get)(v, i))) {
			return 0;
		}
	}
	return 1;
}

VECTOR_T* MATLAB_NAMESPACE::all(const MATRIX_T* m, int dim) {
	if (dim == 1) {
		VECTOR_T* all_v = VECTOR_ID(alloc)(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			VECTOR_ID(const_view) m_col_i = MATRIX_ID(const_column)(m, i);
			VECTOR_ID(set)(all_v, i, all(&m_col_i.vector));
		}
		return all_v;
	} else if (dim == 2) {
		VECTOR_T* all_v = VECTOR_ID(alloc)(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			VECTOR_ID(const_view) m_row_i = MATRIX_ID(const_row)(m, i);
			VECTOR_ID(set)(all_v, i, all(&m_row_i.vector));
		}
		return all_v;
	} else {
		return NULL;
	}
}

int MATLAB_NAMESPACE::any(const VECTOR_T* v) {
	for (int i = 0; i < (int)v->size; i++) {
		if (fp_nonzero(VECTOR_ID(get)(v, i))) {
			return 1;
		}
	}
	return 0;
}

VECTOR_T* MATLAB_NAMESPACE::any(const MATRIX_T* m, int dim) {
	if (dim == 1) {
		VECTOR_T* any_v = VECTOR_ID(alloc)(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			VECTOR_ID(const_view) m_col_i = MATRIX_ID(const_column)(m, i);
			VECTOR_ID(set)(any_v, i, any(&m_col_i.vector));
		}
		return any_v;
	} else if (dim == 2) {
		VECTOR_T* any_v = VECTOR_ID(alloc)(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			VECTOR_ID(const_view) m_row_i = MATRIX_ID(const_row)(m, i);
			VECTOR_ID(set)(any_v, i, any(&m_row_i.vector));
		}
		return any_v;
	} else {
		return NULL;
	}
}

std::string MATLAB_NAMESPACE::dec2bin(int n) {
	return dec2bin(n, 0);
}

std::string MATLAB_NAMESPACE::dec2bin(int n, int len) {
	if (n < 0) {
		return "";
	}
	int binlen = (int)(std::floor(1.0 + std::log(n) / std::log(2)));
	std::string bin((len > binlen) ? len : binlen, '0');
	for (int i = bin.size() - 1; i >= 0; i--) {
		int remainder = n % 2;
		if (remainder) {
			bin[i] = '1';
		}
		n >>= 1;
	}
	return bin;
}

MATRIX_T* MATLAB_NAMESPACE::diag(const VECTOR_T* v, int k) {
	int i0;
	int j0;
	if (k >= 0) { i0 = 0; j0 = k; }
	else { i0 = -k; j0 = 0; }
	int n = (int)v->size + (int)std::abs(k);
	MATRIX_T* diag_m = MATRIX_ID(calloc)(n, n);
	for (int i = 0; i < (int)v->size; i++) {
		MATRIX_ID(set)(diag_m, i0 + i, j0 + i, VECTOR_ID(get)(v, i));
	}
	return diag_m;
}

VECTOR_T* MATLAB_NAMESPACE::diag(const MATRIX_T* m, int k) {
	if (k <= -(int)m->size1 || k >= (int)m->size2) {
		return NULL;
	}
	int i0;
	int j0;
	if (k >= 0) { i0 = 0; j0 = k; }
	else { i0 = -k; j0 = 0; }
	int n_rows = m->size1 - i0;
	int n_cols = m->size2 - j0;
	int n = (n_rows < n_cols) ? n_rows : n_cols;
	VECTOR_T* diag_v = VECTOR_ID(alloc)(n);
	for (int i = 0; i < n; i++) {
		VECTOR_ID(set)(diag_v, i, MATRIX_ID(get)(m, i0 + i, j0 + i));
	}
	return diag_v;
}

MATRIX_T* MATLAB_NAMESPACE::eye(int size) {
	return eye(size, size);
}

MATRIX_T* MATLAB_NAMESPACE::eye(int size1, int size2) {
	MATRIX_T* eye_m = MATRIX_ID(calloc)(size1, size2);
	VECTOR_ID(view) diag_eye_m = MATRIX_ID(diagonal)(eye_m);
	VECTOR_ID(set_all)(&diag_eye_m.vector, 1.0);
	return eye_m;
}

VECTOR_T* MATLAB_NAMESPACE::find(const VECTOR_T* v, int n, const std::string& direction) {
	int n_find = nnz(v);
	if (n_find == 0 || n < 1) {
		return NULL;
	}
	VECTOR_T* find_v = VECTOR_ID(alloc)((n < n_find) ? n : n_find);
	if (direction == "first") {
		int position = 0;
		for (int i = 0; i < (int)v->size && position < (int)find_v->size; i++) {
			if (fp_nonzero(VECTOR_ID(get)(v, i))) {
				VECTOR_ID(set)(find_v, position, i);
				position++;
			}
		}
		return find_v;
	} else if (direction == "last") {
		int position = find_v->size - 1;
		for (int i = v->size - 1; i >= 0 && position >= 0; i--) {
			if (fp_nonzero(VECTOR_ID(get)(v, i))) {
				VECTOR_ID(set)(find_v, position, i);
				position--;
			}
		}
		return find_v;
	} else {
		VECTOR_ID(free)(find_v);
		return NULL;
	}
}

VECTOR_T* MATLAB_NAMESPACE::find(const MATRIX_T* m, int n, const std::string& direction) {
	VECTOR_T* v = to_vector(m);
	VECTOR_T* find_v = find(v, n, direction);
	VECTOR_ID(free)(v);
	return find_v;
}

/*
 * Emulates the two-return version of "find".
 */
MATRIX_T* MATLAB_NAMESPACE::find_ij(const MATRIX_T* m, int n, const std::string& direction) {
	VECTOR_T* find_v = find(m, n, direction);
	if (find_v == NULL) {
		return NULL;
	} else {
		MATRIX_T* find_m = MATRIX_ID(alloc)(find_v->size, 2);
		for (int i = 0; i < (int)find_v->size; i++) {
			int index = (int)VECTOR_ID(get)(find_v, i);
			int row = index % (int)m->size1;
			int column = index / (int)m->size1;
			MATRIX_ID(set)(find_m, i, 0, (FP_T)row);
			MATRIX_ID(set)(find_m, i, 1, (FP_T)column);
		}
		VECTOR_ID(free)(find_v);
		return find_m;
	}
}

VECTOR_T* MATLAB_NAMESPACE::hist(const VECTOR_T* v, int n) {
	VECTOR_T* centers = VECTOR_ID(alloc)(n);
	FP_T min_value = min(v);
	FP_T max_value = max(v);
	FP_T width = (max_value - min_value) / (FP_T)n;
	for (int i = 0; i < n; i++) {
		VECTOR_ID(set)(centers, i, min_value + (i + 0.5) * width);
	}
	VECTOR_T* hist_v = hist(v, centers);
	VECTOR_ID(free)(centers);
	return hist_v;
}

VECTOR_T* MATLAB_NAMESPACE::hist(const VECTOR_T* v, const VECTOR_T* centers) {
	int n = centers->size;
	VECTOR_T* hist_v = VECTOR_ID(calloc)(n);
	for (int i = 0; i < (int)v->size; i++) {
		FP_T value = VECTOR_ID(get)(v, i);
		int index = n - 1;
		for (int j = 0; j < n - 1; j++) {
			FP_T left = VECTOR_ID(get)(centers, j);
			FP_T right = VECTOR_ID(get)(centers, j + 1);
			if (value < left) {
				index = j;
				break;
			} else if (value < right) {
				FP_T middle = (left + right) / 2.0;
				if (fp_less_or_equal(value, middle)) {
					index = j;
				} else {
					index = j + 1;
				}
				break;
			}
		}
		VECTOR_ID(set)(hist_v, index, VECTOR_ID(get)(hist_v, index) + 1.0);
	}
	return hist_v;
}

MATRIX_T* MATLAB_NAMESPACE::inv(const MATRIX_T* m) {
	if (m->size1 != m->size2) {
		return NULL;
	}
	gsl_matrix* lu = to_matrix_double(m);
	gsl_permutation* p = gsl_permutation_alloc(m->size1);
	int signum;
	gsl_linalg_LU_decomp(lu, p, &signum);
	gsl_matrix* inv_m = NULL;
	if (fp_nonzero(gsl_linalg_LU_det(lu, signum))) {
		inv_m = gsl_matrix_alloc(m->size1, m->size2);
		gsl_linalg_LU_invert(lu, p, inv_m);
	}
	gsl_matrix_free(lu);
	gsl_permutation_free(p);
	if (inv_m == NULL) {
		return NULL;
	} else {
		MATRIX_T* ret = to_matrix(inv_m);
		gsl_matrix_free(inv_m);
		return ret;
	}
}

int MATLAB_NAMESPACE::length(const VECTOR_T* v) {
	return v->size;
}

int MATLAB_NAMESPACE::length(const MATRIX_T* m) {
	return (m->size1 > m->size2) ? m->size1 : m->size2;
}

FP_T MATLAB_NAMESPACE::max(FP_T x, FP_T y) {
	if (gsl_isnan((double)x) == 1) {
		return y;
	} else if (gsl_isnan((double)y) == 1) {
		return x;
	} else {
		return (x > y) ? x : y;
	}
}

FP_T MATLAB_NAMESPACE::max(const VECTOR_T* v) {
	FP_T max = (FP_T)GSL_NAN;
	int i = 0;
	for ( ; i < (int)v->size; i++) {
		FP_T value = VECTOR_ID(get)(v, i);
		if (gsl_isnan((double)value) == 0) {
			max = value;
			break;
		}
	}
	for ( ; i < (int)v->size; i++) {
		FP_T value = VECTOR_ID(get)(v, i);
		if (gsl_isnan((double)value) == 0 && value > max) {
			max = value;
		}
	}
	return max;
}

/*
 * Emulates (max(m)) or (max(m')).
 */
VECTOR_T* MATLAB_NAMESPACE::max(const MATRIX_T* m, int dim) {
	if (dim == 1) {
		VECTOR_T* max_v = VECTOR_ID(alloc)(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			VECTOR_ID(const_view) m_col_i = MATRIX_ID(const_column)(m, i);
			FP_T value = max(&m_col_i.vector);
			VECTOR_ID(set)(max_v, i, value);
		}
		return max_v;
	} else if (dim == 2) {
		VECTOR_T* max_v = VECTOR_ID(alloc)(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			VECTOR_ID(const_view) m_row_i = MATRIX_ID(const_row)(m, i);
			FP_T value = max(&m_row_i.vector);
			VECTOR_ID(set)(max_v, i, value);
		}
		return max_v;
	} else {
		return NULL;
	}
}

FP_T MATLAB_NAMESPACE::mean(const VECTOR_T* v, const std::string& opt) {
	if (opt == "a") {
		FP_T sum = 0.0;
		for (int i = 0; i < (int)v->size; i++) {
			sum += VECTOR_ID(get)(v, i);
		}
		return sum / (FP_T)v->size;
	} else if (opt == "g") {
		FP_T product = 1.0;
		for (int i = 0; i < (int)v->size; i++) {
			product *= VECTOR_ID(get)(v, i);
		}
		return std::pow(product, (FP_T)1.0 / (FP_T)v->size);
	} else if (opt == "h") {
		FP_T sum = 0.0;
		for (int i = 0; i < (int)v->size; i++) {
			sum += 1.0 / VECTOR_ID(get)(v, i);
		}
		return (FP_T)v->size / sum;
	} else {
		return GSL_NAN;
	}
}

VECTOR_T* MATLAB_NAMESPACE::mean(const MATRIX_T* m, int dim, const std::string& opt) {
	if (dim == 1) {
		VECTOR_T* mean_v = VECTOR_ID(alloc)(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			VECTOR_ID(const_view) m_col_i = MATRIX_ID(const_column)(m, i);
			FP_T value = mean(&m_col_i.vector, opt);
			VECTOR_ID(set)(mean_v, i, value);
		}
		return mean_v;
	} else if (dim == 2) {
		VECTOR_T* mean_v = VECTOR_ID(alloc)(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			VECTOR_ID(const_view) m_row_i = MATRIX_ID(const_row)(m, i);
			FP_T value = mean(&m_row_i.vector, opt);
			VECTOR_ID(set)(mean_v, i, value);
		}
		return mean_v;
	} else {
		return NULL;
	}
}

FP_T MATLAB_NAMESPACE::min(FP_T x, FP_T y) {
	if (gsl_isnan((double)x) == 1) {
		return y;
	} else if (gsl_isnan((double)y) == 1) {
		return x;
	} else {
		return (x < y) ? x : y;
	}
}

FP_T MATLAB_NAMESPACE::min(const VECTOR_T* v) {
	FP_T min = (FP_T)GSL_NAN;
	int i = 0;
	for ( ; i < (int)v->size; i++) {
		FP_T value = VECTOR_ID(get)(v, i);
		if (gsl_isnan((double)value) == 0) {
			min = value;
			break;
		}
	}
	for ( ; i < (int)v->size; i++) {
		FP_T value = VECTOR_ID(get)(v, i);
		if (gsl_isnan((double)value) == 0 && value < min) {
			min = value;
		}
	}
	return min;
}

/*
 * Emulates (min(m)) or (min(m')).
 */
VECTOR_T* MATLAB_NAMESPACE::min(const MATRIX_T* m, int dim) {
	if (dim == 1) {
		VECTOR_T* min_v = VECTOR_ID(alloc)(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			VECTOR_ID(const_view) m_col_i = MATRIX_ID(const_column)(m, i);
			FP_T value = min(&m_col_i.vector);
			VECTOR_ID(set)(min_v, i, value);
		}
		return min_v;
	} else if (dim == 2) {
		VECTOR_T* min_v = VECTOR_ID(alloc)(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			VECTOR_ID(const_view) m_row_i = MATRIX_ID(const_row)(m, i);
			FP_T value = min(&m_row_i.vector);
			VECTOR_ID(set)(min_v, i, value);
		}
		return min_v;
	} else {
		return NULL;
	}
}

int MATLAB_NAMESPACE::nnz(const VECTOR_T* v) {
	int nnz = 0;
	for (int i = 0; i < (int)v->size; i++) {
		if (fp_nonzero(VECTOR_ID(get)(v, i))) {
			nnz++;
		}
	}
	return nnz;
}

int MATLAB_NAMESPACE::nnz(const MATRIX_T* m) {
	VECTOR_T* v = to_vector(m);
	int nnz_v = nnz(v);
	VECTOR_ID(free)(v);
	return nnz_v;
}

VECTOR_T* MATLAB_NAMESPACE::nonzeros(const MATRIX_T* m) {
	VECTOR_T* nz_v = find(m);
	if (nz_v != NULL) {
		for (int i = 0; i < (int)nz_v->size; i++) {
			int i_m = (int)VECTOR_ID(get)(nz_v, i);
			FP_T value = ordinal_index(m, i_m);
			VECTOR_ID(set)(nz_v, i, value);
		}
	}
	return nz_v;
}

/*
 * Currently only supports p > 1.
 */
FP_T MATLAB_NAMESPACE::norm(const VECTOR_T* v, int p) {
	if (p > 1) {
		FP_T sum = 0.0;
		for (int i = 0; i < (int)v->size; i++) {
			sum += std::pow(std::abs(VECTOR_ID(get)(v, i)), p);
		}
		return std::pow(sum, (FP_T)1.0 / p);
	} else {
		return 0.0;
	}
}

VECTOR_T* MATLAB_NAMESPACE::normpdf(const VECTOR_T* v, FP_T mean, FP_T stdev) {
	VECTOR_T* pdf_v = VECTOR_ID(alloc)(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		double x = (double)VECTOR_ID(get)(v, i);
		double p = gsl_ran_gaussian_pdf(x - mean, stdev);
		VECTOR_ID(set)(pdf_v, i, (FP_T)p);
	}
	return pdf_v;
}

MATRIX_T* MATLAB_NAMESPACE::ones(int size) {
	return ones(size, size);
}

MATRIX_T* MATLAB_NAMESPACE::ones(int size1, int size2) {
	MATRIX_T* ones_m = MATRIX_ID(alloc)(size1, size2);
	MATRIX_ID(set_all)(ones_m, 1.0);
	return ones_m;
}

/*
 * Emulates (ones(size, 1)) or (ones(1, size)).
 */
VECTOR_T* MATLAB_NAMESPACE::ones_vector(int size) {
	VECTOR_T* ones_v = VECTOR_ID(alloc)(size);
	VECTOR_ID(set_all)(ones_v, 1.0);
	return ones_v;
}

FP_T MATLAB_NAMESPACE::prod(const VECTOR_T* v) {
	FP_T prod = 1.0;
	for (int i = 0; i < (int)v->size; i++) {
		prod *= VECTOR_ID(get)(v, i);
	}
	return prod;
}

VECTOR_T* MATLAB_NAMESPACE::prod(const MATRIX_T* m, int dim) {
	if (dim == 1) {
		VECTOR_T* prod_v = VECTOR_ID(alloc)(m->size2);
		VECTOR_ID(set_all)(prod_v, 1.0);
		for (int i = 0; i < (int)m->size1; i++) {
			VECTOR_ID(const_view) m_row_i = MATRIX_ID(const_row)(m, i);
			VECTOR_ID(mul)(prod_v, &m_row_i.vector);
		}
		return prod_v;
	} else if (dim == 2) {
		VECTOR_T* prod_v = VECTOR_ID(alloc)(m->size1);
		VECTOR_ID(set_all)(prod_v, 1.0);
		for (int i = 0; i < (int)m->size2; i++) {
			VECTOR_ID(const_view) m_col_i = MATRIX_ID(const_column)(m, i);
			VECTOR_ID(mul)(prod_v, &m_col_i.vector);
		}
		return prod_v;
	} else {
		return NULL;
	}
}

MATRIX_T* MATLAB_NAMESPACE::rand(int size) {
	return rand(size, size);
}

MATRIX_T* MATLAB_NAMESPACE::rand(int size1, int size2) {
	gsl_rng* rng = get_rng();
	MATRIX_T* rand_m = MATRIX_ID(alloc)(size1, size2);
	for (int i = 0; i < size1; i++) {
		for (int j = 0; j < size2; j++) {
			MATRIX_ID(set)(rand_m, i, j, (FP_T)gsl_rng_uniform(rng));
		}
	}
	return rand_m;
}

VECTOR_T* MATLAB_NAMESPACE::rand_vector(int size) {
	gsl_rng* rng = get_rng();
	VECTOR_T* rand_v = VECTOR_ID(alloc)(size);
	for (int i = 0; i < size; i++) {
		VECTOR_ID(set)(rand_v, i, (FP_T)gsl_rng_uniform(rng));
	}
	return rand_v;
}

/*
 * Generates a permutation of the integers 0 to (size - 1), whereas the MATLAB
 * version uses the integers 1 to size.
 */
gsl_permutation* MATLAB_NAMESPACE::randperm(int size) {
	gsl_rng* rng = get_rng();
	FP_T values[size];
	for (int i = 0; i < size; i++) {
		values[i] = (FP_T)i;
	}
	gsl_ran_shuffle(rng, values, size, sizeof(FP_T));
	VECTOR_ID(view) values_vv = VECTOR_ID(view_array)(values, size);
	gsl_permutation* values_p = to_permutation(&values_vv.vector);
	return values_p;
}

VECTOR_T* MATLAB_NAMESPACE::reverse(const VECTOR_T* v) {
	VECTOR_T* rev_v = VECTOR_ID(alloc)(v->size);
	for (int i = 0; i < (int)v->size; i++) {
		VECTOR_ID(set)(rev_v, i, VECTOR_ID(get)(v, v->size - 1 - i));
	}
	return rev_v;
}

VECTOR_T* MATLAB_NAMESPACE::setxor(const VECTOR_T* v1, const VECTOR_T* v2) {
	if (v1 == NULL && v2 == NULL) {
		return NULL;
	}
	if (v1 == NULL) {
		return unique(v2);
	}
	if (v2 == NULL) {
		return unique(v1);
	}
	VECTOR_T* unique_v1 = unique(v1);
	VECTOR_T* unique_v2 = unique(v2);
	VECTOR_T* unsized_v = VECTOR_ID(alloc)(v1->size + v2->size);
	int n = 0;
	for (int i = 0; i < (int)unique_v1->size; i++) {
		bool found = false;
		FP_T v1_value = VECTOR_ID(get)(unique_v1, i);
		for (int j = 0; j < (int)unique_v2->size; j++) {
			FP_T v2_value = VECTOR_ID(get)(unique_v2, j);
			if (fp_equal(v1_value, v2_value)) {
				found = true;
				break;
			}
		}
		if (!found) {
			VECTOR_ID(set)(unsized_v, n++, v1_value);
		}
	}
	for (int i = 0; i < (int)unique_v2->size; i++) {
		bool found = false;
		FP_T v2_value = VECTOR_ID(get)(unique_v2, i);
		for (int j = 0; j < (int)unique_v1->size; j++) {
			FP_T v1_value = VECTOR_ID(get)(unique_v1, j);
			if (fp_equal(v2_value, v1_value)) {
				found = true;
				break;
			}
		}
		if (!found) {
			VECTOR_ID(set)(unsized_v, n++, v2_value);
		}
	}
	VECTOR_ID(free)(unique_v1);
	VECTOR_ID(free)(unique_v2);
	if (n > 0) {
		VECTOR_T* unsorted_v = VECTOR_ID(alloc)(n);
		VECTOR_ID(view) unsized_subv = VECTOR_ID(subvector)(unsized_v, 0, n);
		VECTOR_ID(memcpy)(unsorted_v, &unsized_subv.vector);
		VECTOR_ID(free)(unsized_v);
		VECTOR_T* setxor_v = sort(unsorted_v);
		VECTOR_ID(free)(unsorted_v);
		return setxor_v;
	} else {
		VECTOR_ID(free)(unsized_v);
		return NULL;
	}
}

VECTOR_T* MATLAB_NAMESPACE::sort(const VECTOR_T* v, const std::string& mode, VECTOR_T** ind) {
	if (mode != "ascend" && mode != "descend") {
		return NULL;
	}
	FP_T elements[v->size];
	to_array(v, elements);
	std::size_t indices[v->size];
	if (mode == "ascend") {
		stable_sort_index(indices, elements, v->size);
	} else {
		stable_sort_index(indices, elements, v->size, fp_greater);
	}
	VECTOR_T* sort_v = VECTOR_ID(alloc)(v->size);
	if (ind != NULL) {
		*ind = VECTOR_ID(alloc)(v->size);
	}
	for (int i = 0; i < (int)v->size; i++) {
		int index = indices[i];
		VECTOR_ID(set)(sort_v, i, elements[index]);
		if (ind != NULL) {
			VECTOR_ID(set)(*ind, i, (FP_T)index);
		}
	}
	return sort_v;
}

MATRIX_T* MATLAB_NAMESPACE::sort(const MATRIX_T* m, int dim, const std::string& mode, MATRIX_T** ind) {
	if (mode != "ascend" && mode != "descend") {
		return NULL;
	}
	if (dim == 1) {
		MATRIX_T* sort_m = MATRIX_ID(alloc)(m->size1, m->size2);
		if (ind != NULL) {
			*ind = MATRIX_ID(alloc)(m->size1, m->size2);
		}
		for (int i = 0; i < (int)m->size2; i++) {
			VECTOR_ID(const_view) m_col_i = MATRIX_ID(const_column)(m, i);
			VECTOR_T* sort_m_col_i;
			if (ind == NULL) {
				sort_m_col_i = sort(&m_col_i.vector, mode);
			} else {
				VECTOR_T* ind_col_i;
				sort_m_col_i = sort(&m_col_i.vector, mode, &ind_col_i);
				MATRIX_ID(set_col)(*ind, i, ind_col_i);
				VECTOR_ID(free)(ind_col_i);
			}
			MATRIX_ID(set_col)(sort_m, i, sort_m_col_i);
			VECTOR_ID(free)(sort_m_col_i);
		}
		return sort_m;
	} else if (dim == 2) {
		MATRIX_T* sort_m = MATRIX_ID(alloc)(m->size1, m->size2);
		if (ind != NULL) {
			*ind = MATRIX_ID(alloc)(m->size1, m->size2);
		}
		for (int i = 0; i < (int)m->size1; i++) {
			VECTOR_ID(const_view) m_row_i = MATRIX_ID(const_row)(m, i);
			VECTOR_T* sort_m_row_i;
			if (ind == NULL) {
				sort_m_row_i = sort(&m_row_i.vector, mode);
			} else {
				VECTOR_T* ind_row_i;
				sort_m_row_i = sort(&m_row_i.vector, mode, &ind_row_i);
				MATRIX_ID(set_row)(*ind, i, ind_row_i);
				VECTOR_ID(free)(ind_row_i);
			}
			MATRIX_ID(set_row)(sort_m, i, sort_m_row_i);
			VECTOR_ID(free)(sort_m_row_i);
		}
		return sort_m;
	} else {
		return NULL;
	}
}

/*
 * Emulates (sortrows(v)) for a column vector.
 */
VECTOR_T* MATLAB_NAMESPACE::sortrows(const VECTOR_T* v, VECTOR_T** ind) {
	return sort(v, "ascend", ind);
}

MATRIX_T* MATLAB_NAMESPACE::sortrows(const MATRIX_T* m, VECTOR_T** ind) {
	VECTOR_T* rows[m->size1];
	for (int i = 0; i < (int)m->size1; i++) {
		rows[i] = VECTOR_ID(alloc)(m->size2);
		MATRIX_ID(get_row)(rows[i], m, i);
	}
	std::size_t indices[m->size1];
	stable_sort_index(indices, rows, m->size1, vector_less);
	for (int i = 0; i < (int)m->size1; i++) {
		VECTOR_ID(free)(rows[i]);
	}
	MATRIX_T* sort_m = MATRIX_ID(alloc)(m->size1, m->size2);
	if (ind != NULL) {
		*ind = VECTOR_ID(alloc)(m->size1);
	}
	for (int i = 0; i < (int)m->size1; i++) {
		int index = indices[i];
		VECTOR_ID(const_view) m_row_index = MATRIX_ID(const_row)(m, index);
		MATRIX_ID(set_row)(sort_m, i, &m_row_index.vector);
		if (ind != NULL) {
			VECTOR_ID(set)(*ind, i, (FP_T)index);
		}
	}
	return sort_m;
}

FP_T MATLAB_NAMESPACE::std(const VECTOR_T* v, int opt) {
	FP_T mu = mean(v);
	FP_T err = 0.0;
	for (int i = 0; i < (int)v->size; i++) {
		err += std::pow(VECTOR_ID(get)(v, i) - mu, 2);
	}
	if (opt == 0) {
		return std::sqrt(err / (FP_T)(v->size - 1));
	} else if (opt == 1) {
		return std::sqrt(err / (FP_T)v->size);
	} else {
		return GSL_NAN;
	}
}

VECTOR_T* MATLAB_NAMESPACE::std(const MATRIX_T* m, int opt, int dim) {
	if (dim == 1) {
		VECTOR_T* std_v = VECTOR_ID(alloc)(m->size2);
		for (int i = 0; i < (int)m->size2; i++) {
			VECTOR_ID(const_view) m_col_i = MATRIX_ID(const_column)(m, i);
			FP_T value = MATLAB_NAMESPACE::std(&m_col_i.vector, opt);
			VECTOR_ID(set)(std_v, i, value);
		}
		return std_v;
	} else if (dim == 2) {
		VECTOR_T* std_v = VECTOR_ID(alloc)(m->size1);
		for (int i = 0; i < (int)m->size1; i++) {
			VECTOR_ID(const_view) m_row_i = MATRIX_ID(const_row)(m, i);
			FP_T value = MATLAB_NAMESPACE::std(&m_row_i.vector, opt);
			VECTOR_ID(set)(std_v, i, value);
		}
		return std_v;
	} else {
		return NULL;
	}
}

FP_T MATLAB_NAMESPACE::sum(const VECTOR_T* v) {
	FP_T sum = 0.0;
	for (int i = 0; i < (int)v->size; i++) {
		sum += VECTOR_ID(get)(v, i);
	}
	return sum;
}

VECTOR_T* MATLAB_NAMESPACE::sum(const MATRIX_T* m, int dim) {
	if (dim == 1) {
		VECTOR_T* sum_v = VECTOR_ID(calloc)(m->size2);
		for (int i = 0; i < (int)m->size1; i++) {
			VECTOR_ID(const_view) m_row_i = MATRIX_ID(const_row)(m, i);
			VECTOR_ID(add)(sum_v, &m_row_i.vector);
		}
		return sum_v;
	} else if (dim == 2) {
		VECTOR_T* sum_v = VECTOR_ID(calloc)(m->size1);
		for (int i = 0; i < (int)m->size2; i++) {
			VECTOR_ID(const_view) m_col_i = MATRIX_ID(const_column)(m, i);
			VECTOR_ID(add)(sum_v, &m_col_i.vector);
		}
		return sum_v;
	} else {
		return NULL;
	}
}

MATRIX_T* MATLAB_NAMESPACE::toeplitz(const VECTOR_T* column, const VECTOR_T* row) {
	const VECTOR_T* _row;
	if (row == NULL) {
		_row = column;
	} else {
		_row = row;
	}
	MATRIX_T* toe_m = MATRIX_ID(alloc)(column->size, _row->size);
	for (int i = 0; i < (int)column->size; i++) {
		for (int j = 0; j < (int)_row->size; j++) {
			FP_T value;
			if (i - j >= 0) {
				value = VECTOR_ID(get)(column, i - j);
			} else {
				value = VECTOR_ID(get)(_row, j - i);
			}
			MATRIX_ID(set)(toe_m, i, j, value);
		}
	}
	return toe_m;
}

MATRIX_T* MATLAB_NAMESPACE::tril(const MATRIX_T* m, int k) {
	if (k < -(int)m->size1 || k > (int)m->size2) {
		return NULL;
	}
	MATRIX_T* tril_m = copy(m);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = i + k + 1; j < (int)m->size2; j++) {
			if (j >= 0) {
				MATRIX_ID(set)(tril_m, i, j, 0.0);
			}
		}
	}
	return tril_m;
}

MATRIX_T* MATLAB_NAMESPACE::triu(const MATRIX_T* m, int k) {
	if (k < -(int)m->size1 || k > (int)m->size2) {
		return NULL;
	}
	MATRIX_T* triu_m = copy(m);
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = i + k - 1; j >= 0; j--) {
			if (j < (int)m->size2) {
				MATRIX_ID(set)(triu_m, i, j, 0.0);
			}
		}
	}
	return triu_m;
}

VECTOR_T* MATLAB_NAMESPACE::unique(const VECTOR_T* v, const std::string& first_or_last, VECTOR_T** i, VECTOR_T** j) {
	if (first_or_last != "first" && first_or_last != "last") {
		return NULL;
	}
	VECTOR_T* sort_v = sort(v);
	VECTOR_T* unsized_v = VECTOR_ID(alloc)(v->size);
	VECTOR_ID(set)(unsized_v, 0, VECTOR_ID(get)(sort_v, 0));
	int n = 1;
	for (int x = 1; x < (int)v->size; x++) {
		FP_T prev_value = VECTOR_ID(get)(sort_v, x - 1);
		FP_T value = VECTOR_ID(get)(sort_v, x);
		if (fp_not_equal(prev_value, value)) {
			VECTOR_ID(set)(unsized_v, n++, value);
		}
	}
	VECTOR_ID(free)(sort_v);
	VECTOR_T* unique_v = VECTOR_ID(alloc)(n);
	VECTOR_ID(view) unsized_subv = VECTOR_ID(subvector)(unsized_v, 0, n);
	VECTOR_ID(memcpy)(unique_v, &unsized_subv.vector);
	VECTOR_ID(free)(unsized_v);
	if (i != NULL) {
		*i = VECTOR_ID(alloc)(n);
		for (int x = 0; x < n; x++) {
			for (int y = 0; y < (int)v->size; y++) {
				if (fp_equal(VECTOR_ID(get)(unique_v, x), VECTOR_ID(get)(v, y))) {
					VECTOR_ID(set)(*i, x, y);
					if (first_or_last == "first") {
						break;
					}
				}
			}
		}
	}
	if (j != NULL) {
		*j = VECTOR_ID(alloc)(v->size);
		for (int x = 0; x < (int)v->size; x++) {
			for (int y = 0; y < n; y++) {
				if (fp_equal(VECTOR_ID(get)(v, x), VECTOR_ID(get)(unique_v, y))) {
					VECTOR_ID(set)(*j, x, y);
					break;
				}
			}
		}
	}
	return unique_v;
}

VECTOR_T* MATLAB_NAMESPACE::unique(const MATRIX_T* m, const std::string& first_or_last, VECTOR_T** i, VECTOR_T** j) {
	VECTOR_T* v = to_vector(m);
	VECTOR_T* unique_v = unique(v, first_or_last, i, j);
	VECTOR_ID(free)(v);
	return unique_v;
}

/*
 * Emulates (unique(m, "rows", first_or_last)).
 */
MATRIX_T* MATLAB_NAMESPACE::unique_rows(const MATRIX_T* m, const std::string& first_or_last, VECTOR_T** i, VECTOR_T** j) {
	if (first_or_last != "first" && first_or_last != "last") {
		return NULL;
	}
	MATRIX_T* sort_m = sortrows(m);
	MATRIX_T* unsized_m = MATRIX_ID(alloc)(m->size1, m->size2);
	VECTOR_ID(view) first_row = MATRIX_ID(row)(sort_m, 0);
	MATRIX_ID(set_row)(unsized_m, 0, &first_row.vector);
	int n_unique = 1;
	for (int x = 1; x < (int)m->size1; x++) {
		VECTOR_ID(view) sort_m_row_x_sub_1 = MATRIX_ID(row)(sort_m, x - 1);
		VECTOR_ID(view) sort_m_row_x = MATRIX_ID(row)(sort_m, x);
		if (compare_vectors(&sort_m_row_x_sub_1.vector, &sort_m_row_x.vector) != 0) {
			MATRIX_ID(set_row)(unsized_m, n_unique++, &sort_m_row_x.vector);
		}
	}
	MATRIX_ID(free)(sort_m);
	MATRIX_T* unique_m = MATRIX_ID(alloc)(n_unique, m->size2);
	MATRIX_ID(view) unsized_subm = MATRIX_ID(submatrix)(unsized_m, 0, 0, n_unique, m->size2);
	MATRIX_ID(memcpy)(unique_m, &unsized_subm.matrix);
	MATRIX_ID(free)(unsized_m);
	if (i != NULL) {
		*i = VECTOR_ID(alloc)(n_unique);
		for (int x = 0; x < n_unique; x++) {
			VECTOR_ID(view) unique_m_row_x = MATRIX_ID(row)(unique_m, x);
			for (int y = 0; y < (int)m->size1; y++) {
				VECTOR_ID(const_view) m_row_y = MATRIX_ID(const_row)(m, y);
				if (compare_vectors(&unique_m_row_x.vector, &m_row_y.vector) == 0) {
					VECTOR_ID(set)(*i, x, y);
					if (first_or_last == "first") {
						break;
					}
				}
			}
		}
	}
	if (j != NULL) {
		*j = VECTOR_ID(alloc)(m->size1);
		for (int x = 0; x < (int)m->size1; x++) {
			VECTOR_ID(const_view) m_row_x = MATRIX_ID(const_row)(m, x);
			for (int y = 0; y < n_unique; y++) {
				VECTOR_ID(view) unique_m_row_y = MATRIX_ID(row)(unique_m, y);
				if (compare_vectors(&m_row_x.vector, &unique_m_row_y.vector) == 0) {
					VECTOR_ID(set)(*j, x, y);
					break;
				}
			}
		}
	}
	return unique_m;
}

MATRIX_T* MATLAB_NAMESPACE::zeros(int size) {
	return MATRIX_ID(calloc)(size, size);
}

MATRIX_T* MATLAB_NAMESPACE::zeros(int size1, int size2) {
	return MATRIX_ID(calloc)(size1, size2);
}

/*
 * Emulates (zeros(size, 1)) or (zeros(1, size)).
 */
VECTOR_T* MATLAB_NAMESPACE::zeros_vector(int size) {
	return VECTOR_ID(calloc)(size);
}
