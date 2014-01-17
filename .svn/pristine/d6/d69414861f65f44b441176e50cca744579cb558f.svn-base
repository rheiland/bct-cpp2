#include <bct/bct_long_double.h>
#include <cstdio>
#include <gsl/gsl_matrix_long_double.h>
#include <gsl/gsl_vector_long_double.h>

int main() {
	
	// Create an uninitialized 30-node square matrix
	gsl_matrix_long_double* m = gsl_matrix_long_double_alloc(30, 30);
	
	// Initialize the matrix with the data in example.dat
	FILE* f = std::fopen("example.dat", "r");
	gsl_matrix_long_double_fscanf(f, m);
	std::fclose(f);
	
	// Display the matrix
	bct_long_double::printf(m, "%Lg");
	
	// Declare variables for optional returns
	gsl_vector_long_double* id;  // In-degree
	gsl_vector_long_double* od;  // Out-degree
	
	// Calculate degree distribution
	gsl_vector_long_double* deg = bct_long_double::degrees_dir(m, &id, &od);
	
	// Display the results
	bct_long_double::printf(id, "%Lg");
	bct_long_double::printf(od, "%Lg");
	bct_long_double::printf(deg, "%Lg");
	
	// Free all memory
	gsl_matrix_long_double_free(m);    // Could use bct::gsl_free(m)
	gsl_vector_long_double_free(id);   // Could use bct::gsl_free(id)
	gsl_vector_long_double_free(od);   // ...
	gsl_vector_long_double_free(deg);
	
	return 0;
}
