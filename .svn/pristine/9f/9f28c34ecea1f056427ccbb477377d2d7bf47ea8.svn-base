#include <gsl/gsl_errno.h>
#include <sstream>

#include "bct.h"

/*
 * Catches GSL errors and throws BCT exceptions.
 */
void BCT_NAMESPACE::gsl_error_handler(const char* reason, const char* file, int line, int gsl_errno) {
	std::stringstream what;
	what << reason << " in " << file << ", line " << line << ".";
	throw bct_exception(what.str());
}

/*
 * Overloaded convenience function for freeing GSL vectors and matrices.
 */
void BCT_NAMESPACE::gsl_free(VECTOR_T* v) { VECTOR_ID(free)(v); }
void BCT_NAMESPACE::gsl_free(MATRIX_T* m) { MATRIX_ID(free)(m); }
void BCT_NAMESPACE::gsl_free(std::vector<MATRIX_T*>& m) {
	for (int i = 0; i < (int)m.size(); i++) {
		if (m[i] != NULL) {
			MATRIX_ID(free)(m[i]);
			m[i] = NULL;
		}
	}
}
void BCT_NAMESPACE::gsl_free(gsl_permutation* p) { gsl_permutation_free(p); }

/*
 * Initializes the BCT library for external use.
 */
void BCT_NAMESPACE::init() {
	gsl_set_error_handler(gsl_error_handler);
}

/*
 * Returns the number of edges in a directed graph.
 */
int BCT_NAMESPACE::number_of_edges_dir(const MATRIX_T* m) {
	return nnz(m);
}

/*
 * Returns the number of edges in an undirected graph.
 */
int BCT_NAMESPACE::number_of_edges_und(const MATRIX_T* m) {
	MATRIX_T* triu_m = triu(m);
	int ret = nnz(triu_m);
	MATRIX_ID(free)(triu_m);
	return ret;
}

/*
 * Returns the number of nodes in a graph.
 */
int BCT_NAMESPACE::number_of_nodes(const MATRIX_T* m) {
	return (int)m->size1;
}
