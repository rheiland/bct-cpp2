#include "bct.h"

/*
 * Returns a randomized graph with equivalent degree sequence to the original
 * weighted directed graph.  On average, each edge is rewired ITER times.  Out-
 * strength is preserved for weighted graphs, while in-strength is not.
 */
MATRIX_T* BCT_NAMESPACE::randmio_dir(const MATRIX_T* R, int ITER) {
	if (safe_mode) check_status(R, SQUARE | DIRECTED, "randmio_dir");
	
	gsl_rng* rng = get_rng();
	
	// [i j]=find(R);
	MATRIX_T* find_R = find_ij(R);
	VECTOR_ID(view) i = MATRIX_ID(column)(find_R, 0);
	VECTOR_ID(view) j = MATRIX_ID(column)(find_R, 1);
	
	// K=length(i);
	int K = length(&i.vector);
	
	// ITER=K*ITER;
	ITER = K * ITER;
	
	MATRIX_T* _R = copy(R);
	
	// for iter=1:ITER
	for (int iter = 1; iter <= ITER; iter++) {
		
		// while 1
		while (true) {
			
			int e1, e2;
			int a, b, c, d;
			
			// while 1
			while (true) {
				
				// e1=ceil(K*rand);
				e1 = gsl_rng_uniform_int(rng, K);
				
				// e2=ceil(K*rand);
				e2 = gsl_rng_uniform_int(rng, K);
				
				// while (e2==e1),
				while (e2 == e1) {
					
					// e2=ceil(K*rand);
					e2 = gsl_rng_uniform_int(rng, K);
				}
				
				// a=i(e1); b=j(e1);
				a = (int)VECTOR_ID(get)(&i.vector, e1);
				b = (int)VECTOR_ID(get)(&j.vector, e1);
				
				// c=i(e2); d=j(e2);
				c = (int)VECTOR_ID(get)(&i.vector, e2);
				d = (int)VECTOR_ID(get)(&j.vector, e2);
				
				// if all(a~=[c d]) && all(b~=[c d]);
				if (a != c && a != d && b != c && b != d) {
					
					// break
					break;
				}
			}
			
			// if ~(R(a,d) || R(c,b))
			if (fp_zero(MATRIX_ID(get)(_R, a, d)) && fp_zero(MATRIX_ID(get)(_R, c, b))) {
				
				// R(a,d)=R(a,b); R(a,b)=0;
				MATRIX_ID(set)(_R, a, d, MATRIX_ID(get)(_R, a, b));
				MATRIX_ID(set)(_R, a, b, 0.0);
				
				// R(c,b)=R(c,d); R(c,d)=0;
				MATRIX_ID(set)(_R, c, b, MATRIX_ID(get)(_R, c, d));
				MATRIX_ID(set)(_R, c, d, 0.0);
				
				// j(e1) = d;
				VECTOR_ID(set)(&j.vector, e1, (FP_T)d);
				
				// j(e2) = b;
				VECTOR_ID(set)(&j.vector, e2, (FP_T)b);
				
				// break;
				break;
			}
		}
	}
	
	MATRIX_ID(free)(find_R);
	return _R;
}
