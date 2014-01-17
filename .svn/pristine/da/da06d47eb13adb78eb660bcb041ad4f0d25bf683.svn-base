#include <cstdio>

#include "bct.h"

/*
 * Prints a vector using the given format for each element.  This is only
 * provided for debugging purposes.  In other cases, use gsl_vector_fprintf.
 */
void BCT_NAMESPACE::printf(const VECTOR_T* v, const std::string& format) {
	for (int i = 0; i < (int)v->size; i++) {
		std::printf(format.c_str(), VECTOR_ID(get)(v, i));
		std::printf(" ");
	}
	std::printf("\n");
}

/*
 * Prints a matrix using the given format for each element.  This is only
 * provided for debugging purposes.  In other cases, use gsl_matrix_fprintf.
 */
void BCT_NAMESPACE::printf(const MATRIX_T* m, const std::string& format) {
	for (int i = 0; i < (int)m->size1; i++) {
		for (int j = 0; j < (int)m->size2; j++) {
			std::printf(format.c_str(), MATRIX_ID(get)(m, i, j));
			std::printf(" ");
		}
		std::printf("\n");
	}
}

/*
 * Prints a permutation using the given format for each element.  This is only
 * provided for debugging purposes.  In other cases, use
 * gsl_permutation_fprintf.
 */
void BCT_NAMESPACE::printf(const gsl_permutation* p, const std::string& format) {
	for (int i = 0; i < (int)p->size; i++) {
		std::printf(format.c_str(), gsl_permutation_get(p, i));
		std::printf(" ");
	}
	std::printf("\n");
}
