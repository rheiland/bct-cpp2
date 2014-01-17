#include "../precision.h"

#undef SKIP

#ifdef GSL_FLOAT
#ifdef MATLAB_FLOAT_H
#define SKIP
#else
#define MATLAB_FLOAT_H
#endif
#endif

#ifdef GSL_DOUBLE
#ifdef MATLAB_H
#define SKIP
#else
#define MATLAB_H
#endif
#endif

#ifdef GSL_LONG_DOUBLE
#ifdef MATLAB_LONG_DOUBLE_H
#define SKIP
#else
#define MATLAB_LONG_DOUBLE_H
#endif
#endif

#ifndef SKIP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <limits>
#include <string>

namespace MATLAB_NAMESPACE {
	
	// Functions
	VECTOR_T* abs(const VECTOR_T* v);
	MATRIX_T* abs(const MATRIX_T* m);
	int all(const VECTOR_T* v);
	VECTOR_T* all(const MATRIX_T* m, int dim = 1);
	int any(const VECTOR_T* v);
	VECTOR_T* any(const MATRIX_T* m, int dim = 1);
	std::string dec2bin(int n);
	std::string dec2bin(int n, int len);
	MATRIX_T* diag(const VECTOR_T* v, int k = 0);
	VECTOR_T* diag(const MATRIX_T* m, int k = 0);
	MATRIX_T* eye(int size);
	MATRIX_T* eye(int size1, int size2);
	VECTOR_T* find(const VECTOR_T* v, int n = std::numeric_limits<int>::max(), const std::string& direction = "first");
	VECTOR_T* find(const MATRIX_T* m, int n = std::numeric_limits<int>::max(), const std::string& direction = "first");
	MATRIX_T* find_ij(const MATRIX_T* m, int n = std::numeric_limits<int>::max(), const std::string& direction = "first");
	VECTOR_T* hist(const VECTOR_T* v, int n = 10);
	VECTOR_T* hist(const VECTOR_T* v, const VECTOR_T* centers);
	MATRIX_T* inv(const MATRIX_T* m);
	int length(const VECTOR_T* v);
	int length(const MATRIX_T* m);
	FP_T max(FP_T x, FP_T y);
	FP_T max(const VECTOR_T* v);
	VECTOR_T* max(const MATRIX_T* m, int dim = 1);
	FP_T mean(const VECTOR_T* v, const std::string& opt = "a");
	VECTOR_T* mean(const MATRIX_T* m, int dim = 1, const std::string& opt = "a");
	FP_T min(FP_T x, FP_T y);
	FP_T min(const VECTOR_T* v);
	VECTOR_T* min(const MATRIX_T* m, int dim = 1);
	int nnz(const VECTOR_T* v);
	int nnz(const MATRIX_T* m);
	VECTOR_T* nonzeros(const MATRIX_T* m);
	FP_T norm(const VECTOR_T* v, int p);
	VECTOR_T* normpdf(const VECTOR_T* v, FP_T mean, FP_T stdev);
	MATRIX_T* ones(int size);
	MATRIX_T* ones(int size1, int size2);
	VECTOR_T* ones_vector(int size);
	FP_T prod(const VECTOR_T* v);
	VECTOR_T* prod(const MATRIX_T* m, int dim = 1);
	MATRIX_T* rand(int size);
	MATRIX_T* rand(int size1, int size2);
	VECTOR_T* rand_vector(int size);
	gsl_permutation* randperm(int size);
	VECTOR_T* reverse(const VECTOR_T* v);
	VECTOR_T* setxor(const VECTOR_T* v1, const VECTOR_T* v2);
	VECTOR_T* sort(const VECTOR_T* v, const std::string& mode = "ascend", VECTOR_T** ind = NULL);
	MATRIX_T* sort(const MATRIX_T* m, int dim = 1, const std::string& mode = "ascend", MATRIX_T** ind = NULL);
	VECTOR_T* sortrows(const VECTOR_T* v, VECTOR_T** ind = NULL);
	MATRIX_T* sortrows(const MATRIX_T* m, VECTOR_T** ind = NULL);
	FP_T std(const VECTOR_T* v, int opt = 0);
	VECTOR_T* std(const MATRIX_T* m, int opt = 0, int dim = 1);
	FP_T sum(const VECTOR_T* v);
	VECTOR_T* sum(const MATRIX_T* m, int dim = 1);
	MATRIX_T* toeplitz(const VECTOR_T* column, const VECTOR_T* row = NULL);
	MATRIX_T* tril(const MATRIX_T* m, int k = 0);
	MATRIX_T* triu(const MATRIX_T* m, int k = 0);
	VECTOR_T* unique(const VECTOR_T* v, const std::string& first_or_last = "last", VECTOR_T** i = NULL, VECTOR_T** j = NULL);
	VECTOR_T* unique(const MATRIX_T* m, const std::string& first_or_last = "last", VECTOR_T** i = NULL, VECTOR_T** j = NULL);
	MATRIX_T* unique_rows(const MATRIX_T* m, const std::string& first_or_last = "last", VECTOR_T** i = NULL, VECTOR_T** j = NULL);
	MATRIX_T* zeros(int size);
	MATRIX_T* zeros(int size1, int size2);
	VECTOR_T* zeros_vector(int size);
	
	// Operators
	VECTOR_T* concatenate(const VECTOR_T* v, FP_T x);
	VECTOR_T* concatenate(FP_T x, const VECTOR_T* v);
	VECTOR_T* concatenate(const VECTOR_T* v1, const VECTOR_T* v2);
	MATRIX_T* concatenate_columns(const VECTOR_T* v1, const VECTOR_T* v2);
	MATRIX_T* concatenate_columns(const MATRIX_T* m, const VECTOR_T* v);
	MATRIX_T* concatenate_columns(const VECTOR_T* v, const MATRIX_T* m);
	MATRIX_T* concatenate_columns(const MATRIX_T* m1, const MATRIX_T* m2);
	MATRIX_T* concatenate_rows(const VECTOR_T* v1, const VECTOR_T* v2);
	MATRIX_T* concatenate_rows(const MATRIX_T* m, const VECTOR_T* v);
	MATRIX_T* concatenate_rows(const VECTOR_T* v, const MATRIX_T* m);
	MATRIX_T* concatenate_rows(const MATRIX_T* m1, const MATRIX_T* m2);
	VECTOR_T* copy(const VECTOR_T* v);
	MATRIX_T* copy(const MATRIX_T* m);
	MATRIX_T* div_left(const MATRIX_T* m1, const MATRIX_T* m2);
	MATRIX_T* div_right(const MATRIX_T* m1, const MATRIX_T* m2);
	VECTOR_T* logical_and(const VECTOR_T* v1, const VECTOR_T* v2);
	MATRIX_T* logical_and(const MATRIX_T* m1, const MATRIX_T* m2);
	VECTOR_T* logical_not(const VECTOR_T* v);
	MATRIX_T* logical_not(const MATRIX_T* m);
	VECTOR_T* logical_or(const VECTOR_T* v1, const VECTOR_T* v2);
	MATRIX_T* logical_or(const MATRIX_T* m1, const MATRIX_T* m2);
	MATRIX_T* mul(const MATRIX_T* m1, const MATRIX_T* m2);
	MATRIX_T* pow(const MATRIX_T* m, int power);
	VECTOR_T* pow_elements(const VECTOR_T* v, FP_T power);
	VECTOR_T* pow_elements(const VECTOR_T* v, const VECTOR_T* powers);
	MATRIX_T* pow_elements(const MATRIX_T* m, FP_T power);
	MATRIX_T* pow_elements(const MATRIX_T* m, const MATRIX_T* powers);
	VECTOR_T* sequence(int start, int end);
	VECTOR_T* sequence(int start, int step, int end);
	
	// Floating-point comparison
	extern FP_T epsilon;
	int fp_compare(FP_T x, FP_T y);
	bool fp_zero(FP_T x);
	bool fp_nonzero(FP_T x);
	bool fp_equal(FP_T x, FP_T y);
	bool fp_not_equal(FP_T x, FP_T y);
	bool fp_less(FP_T x, FP_T y);
	bool fp_less_or_equal(FP_T x, FP_T y);
	bool fp_greater(FP_T x, FP_T y);
	bool fp_greater_or_equal(FP_T x, FP_T y);
	
	// Vector/matrix comparison
	typedef bool (*comparator)(FP_T, FP_T);
	int compare_vectors(const VECTOR_T* v1, const VECTOR_T* v2);
	bool vector_less(VECTOR_T* v1, VECTOR_T* v2);
	int compare_matrices(const MATRIX_T* m1, const MATRIX_T* m2);
	bool matrix_less(MATRIX_T* m1, MATRIX_T* m2);
	VECTOR_T* compare_elements(const VECTOR_T* v, comparator compare, FP_T x);
	VECTOR_T* compare_elements(const VECTOR_T* v1, comparator compare, const VECTOR_T* v2);
	MATRIX_T* compare_elements(const MATRIX_T* m, comparator compare, FP_T x);
	MATRIX_T* compare_elements(const MATRIX_T* m1, comparator compare, const MATRIX_T* m2);
	
	// Vector-by-vector indexing
	VECTOR_T* ordinal_index(const VECTOR_T* v, const VECTOR_T* indices);
	void ordinal_index_assign(VECTOR_T* v, const VECTOR_T* indices, FP_T value);
	void ordinal_index_assign(VECTOR_T* v, const VECTOR_T* indices, const VECTOR_T* values);
	VECTOR_T* logical_index(const VECTOR_T* v, const VECTOR_T* logical_v);
	void logical_index_assign(VECTOR_T* v, const VECTOR_T* logical_v, FP_T value);
	void logical_index_assign(VECTOR_T* v, const VECTOR_T* logical_v, const VECTOR_T* values);
	
	// Matrix-by-integer indexing
	FP_T ordinal_index(const MATRIX_T* m, int index);
	void ordinal_index_assign(MATRIX_T* m, int index, FP_T value);
	
	// Matrix-by-vector indexing
	VECTOR_T* ordinal_index(const MATRIX_T* m, const VECTOR_T* indices);
	void ordinal_index_assign(MATRIX_T* m, const VECTOR_T* indices, FP_T value);
	void ordinal_index_assign(MATRIX_T* m, const VECTOR_T* indices, const VECTOR_T* values);
	VECTOR_T* logical_index(const MATRIX_T* m, const VECTOR_T* logical_v);
	void logical_index_assign(MATRIX_T* m, const VECTOR_T* logical_v, FP_T value);
	void logical_index_assign(MATRIX_T* m, const VECTOR_T* logical_v, const VECTOR_T* values);
	
	// Matrix-by-two-vectors indexing (non-mixed)
	MATRIX_T* ordinal_index(const MATRIX_T* m, const VECTOR_T* rows, const VECTOR_T* columns);
	void ordinal_index_assign(MATRIX_T* m, const VECTOR_T* rows, const VECTOR_T* columns, FP_T value);
	void ordinal_index_assign(MATRIX_T* m, const VECTOR_T* rows, const VECTOR_T* columns, const MATRIX_T* values);
	MATRIX_T* logical_index(const MATRIX_T* m, const VECTOR_T* logical_rows, const VECTOR_T* logical_columns);
	void logical_index_assign(MATRIX_T* m, const VECTOR_T* logical_rows, const VECTOR_T* logical_columns, FP_T value);
	void logical_index_assign(MATRIX_T* m, const VECTOR_T* logical_rows, const VECTOR_T* logical_columns, const MATRIX_T* values);
	
	// Matrix-by-two-vectors indexing (mixed)
	MATRIX_T* ord_log_index(const MATRIX_T* m, const VECTOR_T* rows, const VECTOR_T* logical_columns);
	void ord_log_index_assign(MATRIX_T* m, const VECTOR_T* rows, const VECTOR_T* logical_columns, FP_T value);
	void ord_log_index_assign(MATRIX_T* m, const VECTOR_T* rows, const VECTOR_T* logical_columns, const MATRIX_T* values);
	MATRIX_T* log_ord_index(const MATRIX_T* m, const VECTOR_T* logical_rows, const VECTOR_T* columns);
	void log_ord_index_assign(MATRIX_T* m, const VECTOR_T* logical_rows, const VECTOR_T* columns, FP_T value);
	void log_ord_index_assign(MATRIX_T* m, const VECTOR_T* logical_rows, const VECTOR_T* columns, const MATRIX_T* values);
	
	// Matrix-by-matrix indexing
	MATRIX_T* ordinal_index(const MATRIX_T* m, const MATRIX_T* indices);
	void ordinal_index_assign(MATRIX_T* m, const MATRIX_T* indices, FP_T value);
	void ordinal_index_assign(MATRIX_T* m, const MATRIX_T* indices, const MATRIX_T* values);
	VECTOR_T* logical_index(const MATRIX_T* m, const MATRIX_T* logical_m);
	void logical_index_assign(MATRIX_T* m, const MATRIX_T* logical_m, FP_T value);
	void logical_index_assign(MATRIX_T* m, const MATRIX_T* logical_m, const VECTOR_T* values);
	
	// Vector/matrix conversion
	void to_array(const VECTOR_T* v, FP_T* array);
	bool to_bool(const VECTOR_T* v);
	bool to_bool(const MATRIX_T* m);
	VECTOR_T* to_vector(const gsl_vector* v_d);
	gsl_vector* to_vector_double(const VECTOR_T* v);
	VECTOR_T* to_vector(const MATRIX_T* m);
	MATRIX_T* to_column_matrix(const VECTOR_T* v);
	MATRIX_T* to_row_matrix(const VECTOR_T* v);
	MATRIX_T* to_matrix(const gsl_matrix* m_d);
	gsl_matrix* to_matrix_double(const MATRIX_T* m);
	VECTOR_T* to_vector(const gsl_permutation* p);
	gsl_permutation* to_permutation(const VECTOR_T* v);
	
	// Utility
	gsl_rng* get_rng();
	void seed_rng(const gsl_rng* rng, unsigned long seed);
	VECTOR_T* permute(const gsl_permutation* p, const VECTOR_T* v);
	MATRIX_T* permute_columns(const gsl_permutation* p, const MATRIX_T* m);
	MATRIX_T* permute_rows(const gsl_permutation* p, const MATRIX_T* m);
}

#endif
