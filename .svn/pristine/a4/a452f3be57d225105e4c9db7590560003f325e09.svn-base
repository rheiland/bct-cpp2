#include <bct/bct.h>
#include <cstdio>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

int main() {
	
	// Create an uninitialized 30-node square matrix
	gsl_matrix* m = gsl_matrix_alloc(30, 30);
	
	// Initialize the matrix with the data in example.dat
	FILE* f = std::fopen("example.dat", "r");
	gsl_matrix_fscanf(f, m);
	std::fclose(f);
	
	// Display the matrix
	bct::printf(m, "%g");
	
	// Declare variables for optional returns
	gsl_vector* id;  // In-degree
	gsl_vector* od;  // Out-degree
	
	// Calculate degree distribution
	gsl_vector* deg = bct::degrees_dir(m, &id, &od);
	
	// Display the results
	bct::printf(id, "%g");
	bct::printf(od, "%g");
	bct::printf(deg, "%g");
	
	// Free all memory
	gsl_matrix_free(m);    // Could use bct::gsl_free(m)
	gsl_vector_free(id);   // Could use bct::gsl_free(id)
	gsl_vector_free(od);   // ...
	gsl_vector_free(deg);
	
	return 0;
}
