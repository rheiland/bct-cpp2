#undef FP_T
#undef FP_ID
#undef VECTOR_T
#undef VECTOR_ID
#undef MATRIX_T
#undef MATRIX_ID

#ifdef GSL_FLOAT
#define MATLAB_NAMESPACE matlab_float
#define BCT_NAMESPACE bct_float
#define FP_T float
#define FP_ID(id) id##_##float
#define VECTOR_T gsl_vector_float
#define VECTOR_ID(id) gsl_vector_float##_##id
#define MATRIX_T gsl_matrix_float
#define MATRIX_ID(id) gsl_matrix_float##_##id
#include <gsl/gsl_matrix_float.h>
#include <gsl/gsl_vector_float.h>
#endif

#ifdef GSL_DOUBLE
#define MATLAB_NAMESPACE matlab
#define BCT_NAMESPACE bct
#define FP_T double
#define FP_ID(id) id##_##double
#define VECTOR_T gsl_vector
#define VECTOR_ID(id) gsl_vector##_##id
#define MATRIX_T gsl_matrix
#define MATRIX_ID(id) gsl_matrix##_##id
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#endif

#ifdef GSL_LONG_DOUBLE
#define MATLAB_NAMESPACE matlab_long_double
#define BCT_NAMESPACE bct_long_double
#define FP_T long double
#define FP_ID(id) id##_##long_double
#define VECTOR_T gsl_vector_long_double
#define VECTOR_ID(id) gsl_vector_long_double##_##id
#define MATRIX_T gsl_matrix_long_double
#define MATRIX_ID(id) gsl_matrix_long_double##_##id
#include <gsl/gsl_matrix_long_double.h>
#include <gsl/gsl_vector_long_double.h>
#endif
