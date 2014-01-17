#include "precision.h"

#undef SKIP

#ifdef GSL_FLOAT
#ifdef BCT_FLOAT_H
#define SKIP
#else
#define BCT_FLOAT_H
#endif
#endif

#ifdef GSL_DOUBLE
#ifdef BCT_H
#define SKIP
#else
#define BCT_H
#endif
#endif

#ifdef GSL_LONG_DOUBLE
#ifdef BCT_LONG_DOUBLE_H
#define SKIP
#else
#define BCT_LONG_DOUBLE_H
#endif
#endif

#ifndef SKIP

#include <stdexcept>
#include <vector>

#include "matlab/matlab.h"

namespace BCT_NAMESPACE {
	using namespace MATLAB_NAMESPACE;
	
	class bct_exception : public std::runtime_error {
	public:
		bct_exception(const std::string& what_arg) : std::runtime_error(what_arg) { }
	};

	// Density, degree, and assortativity
	FP_T assortativity_dir(const MATRIX_T* CIJ);
	FP_T assortativity_und(const MATRIX_T* CIJ);
	VECTOR_T* degrees_dir(const MATRIX_T* CIJ, VECTOR_T** id = NULL, VECTOR_T** od = NULL);
	VECTOR_T* degrees_und(const MATRIX_T* CIJ);
	FP_T density_dir(const MATRIX_T* CIJ);
	FP_T density_und(const MATRIX_T* CIJ);
	MATRIX_T* jdegree(const MATRIX_T* CIJ);
	int jdegree_bl(const MATRIX_T* J);
	int jdegree_id(const MATRIX_T* J);
	int jdegree_od(const MATRIX_T* J);
	MATRIX_T* matching_ind(const MATRIX_T* CIJ);
	MATRIX_T* matching_ind_in(const MATRIX_T* CIJ);
	MATRIX_T* matching_ind_out(const MATRIX_T* CIJ);
	VECTOR_T* strengths_dir(const MATRIX_T* CIJ, VECTOR_T** is = NULL, VECTOR_T** os = NULL);
	VECTOR_T* strengths_und(const MATRIX_T* CIJ);

	// Clustering
	VECTOR_T* clustering_coef_bd(const MATRIX_T* A);
	VECTOR_T* clustering_coef_bu(const MATRIX_T* G);
	VECTOR_T* clustering_coef_wd(const MATRIX_T* W);
	VECTOR_T* clustering_coef_wu(const MATRIX_T* W);
	VECTOR_T* efficiency_local(const MATRIX_T* G);

	// Paths, distances, and cycles
	VECTOR_T* breadth(const MATRIX_T* CIJ, int source, VECTOR_T** branch = NULL);
	MATRIX_T* breadthdist(const MATRIX_T* CIJ, MATRIX_T** D = NULL);
	VECTOR_T* charpath_ecc(const MATRIX_T* D, FP_T* radius = NULL, FP_T* diameter = NULL);
	FP_T charpath_lambda(const MATRIX_T* D);
	FP_T capped_charpath_lambda(const MATRIX_T* G);
	FP_T connectivity_length(const MATRIX_T* D);
	VECTOR_T* cycprob_fcyc(const std::vector<MATRIX_T*>& Pq);
	VECTOR_T* cycprob_pcyc(const std::vector<MATRIX_T*>& Pq);
	MATRIX_T* distance_bin(const MATRIX_T* G);
	MATRIX_T* distance_wei(const MATRIX_T* G);
	FP_T efficiency_global(const MATRIX_T* G, const MATRIX_T* D = NULL);
	std::vector<MATRIX_T*> findpaths(const MATRIX_T* CIJ, const VECTOR_T* sources, int qmax, VECTOR_T** plq = NULL, int* qstop = NULL, MATRIX_T** allpths = NULL, MATRIX_T** util = NULL);
	std::vector<MATRIX_T*> findwalks(const MATRIX_T* CIJ, VECTOR_T** wlq = NULL);
	FP_T normalized_path_length(const MATRIX_T* D, FP_T wmax = 1.0);
	FP_T normalized_path_length_m(const MATRIX_T* G, FP_T wmax = 1.0);
	MATRIX_T* reachdist(const MATRIX_T* CIJ, MATRIX_T** D = NULL);

	// Centrality
	VECTOR_T* betweenness_bin(const MATRIX_T* G);
	VECTOR_T* betweenness_wei(const MATRIX_T* G);
	MATRIX_T* edge_betweenness_bin(const MATRIX_T* G, VECTOR_T** BC = NULL);
	MATRIX_T* edge_betweenness_wei(const MATRIX_T* G, VECTOR_T** BC = NULL);
	MATRIX_T* erange(const MATRIX_T* CIJ, FP_T* eta = NULL, MATRIX_T** Eshort = NULL, FP_T* fs = NULL);
	VECTOR_T* eigenvector_centrality(const MATRIX_T* G);

	// Motifs
	enum motif_mode_enum { MILO, SPORNS };
	extern motif_mode_enum motif_mode;
	motif_mode_enum get_motif_mode();
	void set_motif_mode(motif_mode_enum motif_mode);
	std::vector<MATRIX_T*> find_motif34(int m, int n);
	int find_motif34(const MATRIX_T* m);
	VECTOR_T* motif3funct_bin(const MATRIX_T* W, MATRIX_T** F = NULL);
	MATRIX_T* motif3funct_wei(const MATRIX_T* W, MATRIX_T** Q = NULL, MATRIX_T** F = NULL);
	VECTOR_T* motif3funct_wei_v(const MATRIX_T* W, VECTOR_T** Q = NULL, VECTOR_T** F = NULL);
	MATRIX_T* motif3generate(VECTOR_T** ID = NULL, VECTOR_T** N = NULL);
	VECTOR_T* motif3struct_bin(const MATRIX_T* A, MATRIX_T** F = NULL);
	MATRIX_T* motif3struct_wei(const MATRIX_T* W, MATRIX_T** Q = NULL, MATRIX_T** F = NULL);
	VECTOR_T* motif3struct_wei_v(const MATRIX_T* W, VECTOR_T** Q = NULL, VECTOR_T** F = NULL);
	MATRIX_T* motif4generate(VECTOR_T** ID = NULL, VECTOR_T** N = NULL);
	VECTOR_T* motif4funct_bin(const MATRIX_T* W, MATRIX_T** F = NULL);
	MATRIX_T* motif4funct_wei(const MATRIX_T* W, MATRIX_T** Q = NULL, MATRIX_T** F = NULL);
	VECTOR_T* motif4funct_wei_v(const MATRIX_T* W, VECTOR_T** Q = NULL, VECTOR_T** F = NULL);
	VECTOR_T* motif4struct_bin(const MATRIX_T* A, MATRIX_T** F = NULL);
	MATRIX_T* motif4struct_wei(const MATRIX_T* W, MATRIX_T** Q = NULL, MATRIX_T** F = NULL);
	VECTOR_T* motif4struct_wei_v(const MATRIX_T* W, VECTOR_T** Q = NULL, VECTOR_T** F = NULL);

	// Modularity and community structure
	FP_T modularity_dir(const MATRIX_T* A, VECTOR_T** Ci = NULL);
	FP_T modularity_und(const MATRIX_T* A, VECTOR_T** Ci = NULL);
	FP_T modularity_louvain_und(const MATRIX_T* W, VECTOR_T** Ci = NULL, int N = 100);
	VECTOR_T* module_degree_zscore(const MATRIX_T* A, const VECTOR_T* Ci);
	VECTOR_T* participation_coef(const MATRIX_T* A, const VECTOR_T* Ci);
	
	// Synthetic connection networks
	MATRIX_T* makeevenCIJ(int N, int K, int sz_cl);
	MATRIX_T* makefractalCIJ(int mx_lvl, FP_T E, int sz_cl, int* K = NULL);
	MATRIX_T* makelatticeCIJ(int N, int K);
	MATRIX_T* makerandCIJ_bd(int N, int K);
	MATRIX_T* makerandCIJ_bu(int N, int K);
	MATRIX_T* makerandCIJ_wd(int N, int K, FP_T wmin, FP_T wmax);
	MATRIX_T* makerandCIJ_wd_wp(const MATRIX_T* m);
	MATRIX_T* makerandCIJ_wu(int N, int K, FP_T wmin, FP_T wmax);
	MATRIX_T* makerandCIJ_wu_wp(const MATRIX_T* m);
	MATRIX_T* makerandCIJdegreesfixed(const VECTOR_T* in, const VECTOR_T* out);
	MATRIX_T* makerandCIJdegreesfixed(const MATRIX_T* m);
	MATRIX_T* makeringlatticeCIJ(int N, int K);
	MATRIX_T* maketoeplitzCIJ(int N, int K, FP_T s);
	
	// Graph randomization
	MATRIX_T* latmio_dir(const MATRIX_T* R, int ITER);
	MATRIX_T* latmio_dir_connected(const MATRIX_T* R, int ITER);	
	MATRIX_T* latmio_und(const MATRIX_T* R, int ITER);
	MATRIX_T* latmio_und_connected(const MATRIX_T* R, int ITER);
	MATRIX_T* randmio_dir(const MATRIX_T* R, int ITER);
	MATRIX_T* randmio_dir_connected(const MATRIX_T* R, int ITER);
	MATRIX_T* randmio_und(const MATRIX_T* R, int ITER);
	MATRIX_T* randmio_und_connected(const MATRIX_T* R, int ITER);

	// Data sets
	MATRIX_T* get_cat_all();
	MATRIX_T* get_cat_ctx();
	MATRIX_T* get_fve30();
	MATRIX_T* get_fve32();
	MATRIX_T* get_macaque47();
	MATRIX_T* get_macaque71();
	
	// Matrix status checking
	enum status {
		SQUARE = 1, RECTANGULAR = 2,
		UNDIRECTED = 4, DIRECTED = 8,
		BINARY = 16, WEIGHTED = 32,
		POSITIVE = 64, SIGNED = 128,
		NO_LOOPS = 256, LOOPS = 512
	};
	extern bool safe_mode;
	bool get_safe_mode();
	void set_safe_mode(bool safe_mode);
	bool check_status(const MATRIX_T* m, int flags, const std::string& text);
	bool is_square(const MATRIX_T* m);
	bool is_rectangular(const MATRIX_T* m);
	bool is_undirected(const MATRIX_T* m);
	bool is_directed(const MATRIX_T* m);
	bool is_binary(const MATRIX_T* m);
	bool is_weighted(const MATRIX_T* m);
	bool is_positive(const MATRIX_T* m);
	bool is_signed(const MATRIX_T* m);
	bool has_loops(const MATRIX_T* m);
	bool has_no_loops(const MATRIX_T* m);

	// Matrix conversion
	MATRIX_T* invert_elements(const MATRIX_T* m);
	MATRIX_T* remove_loops(const MATRIX_T* m);
	MATRIX_T* to_binary(const MATRIX_T* m);
	MATRIX_T* to_positive(const MATRIX_T* m);
	MATRIX_T* to_undirected_bin(const MATRIX_T* m);
	MATRIX_T* to_undirected_wei(const MATRIX_T* m);
	
	// Utility
	void gsl_error_handler(const char* reason, const char* file, int line, int gsl_errno);
	void gsl_free(VECTOR_T* v);
	void gsl_free(MATRIX_T* m);
	void gsl_free(std::vector<MATRIX_T*>& m);
	void gsl_free(gsl_permutation* p);
	void init();
	int number_of_edges_dir(const MATRIX_T* m);
	int number_of_edges_und(const MATRIX_T* m);
	int number_of_nodes(const MATRIX_T* m);
	MATRIX_T* threshold_absolute(const MATRIX_T* W, FP_T thr);
	MATRIX_T* threshold_proportional_dir(const MATRIX_T* W, FP_T p);
	MATRIX_T* threshold_proportional_und(const MATRIX_T* W, FP_T p);
	
	// Debugging
	void printf(const VECTOR_T* v, const std::string& format);
	void printf(const MATRIX_T* m, const std::string& format);
	void printf(const gsl_permutation* p, const std::string& format);
}

#endif
