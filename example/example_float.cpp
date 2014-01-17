#include <bct/bct_float.h>
#include <cstdio>
#include <gsl/gsl_matrix_float.h>
#include <gsl/gsl_vector_float.h>

int main() {
	
	// Create an uninitialized 30-node square matrix
	gsl_matrix_float* m = gsl_matrix_float_alloc(30, 30);
	
	// Initialize the matrix with the data in example.dat
	FILE* f = std::fopen("example.dat", "r");
	gsl_matrix_float_fscanf(f, m);
	std::fclose(f);
	
	// Display the matrix
	bct_float::printf(m, "%g");
	
	// Declare variables for optional returns
	gsl_vector_float* id;  // In-degree
	gsl_vector_float* od;  // Out-degree
	
	// Calculate degree distribution
	gsl_vector_float* deg = bct_float::degrees_dir(m, &id, &od);
	
	// Display the results
	bct_float::printf(id, "%g");
	bct_float::printf(od, "%g");
	bct_float::printf(deg, "%g");
	
	// Free all memory
	gsl_matrix_float_free(m);    // Could use bct::gsl_free(m)
	gsl_vector_float_free(id);   // Could use bct::gsl_free(id)
	gsl_vector_float_free(od);   // ...
	gsl_vector_float_free(deg);
	
	return 0;
}
