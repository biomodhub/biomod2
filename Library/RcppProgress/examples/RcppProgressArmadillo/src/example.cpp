#include <RcppArmadillo.h>
#include "progress.hpp"


using namespace Rcpp;

// your function for which to provide support
double your_long_computation(int nb) {
	double sum = 0;
	for (int i = 0; i < nb; ++i) {
		if ( Progress::check_abort() )
			return -1;
		for (int j = 0; j < nb; ++j) {
		    arma::mat m1 = arma::eye<arma::mat>(3, 3);
		    arma::mat m2 = arma::eye<arma::mat>(3, 3);

		    m1 = m1 + 3 * (m1 + m2);
			sum += m1[0];
		}
	}

	return sum;
}

void test_sequential2(int max, int nb, bool display_progress) {

	Progress p(max, display_progress);
	for (int i = 0; i < max; ++i) {
		if ( p.increment() ) {
			your_long_computation(nb);
		}
	}
}

void test_multithreaded_omp2(int max, int nb, int threads, bool display_progress) {

#ifdef _OPENMP
	if ( threads > 0 )
		omp_set_num_threads( threads );
	REprintf("Number of threads=%i\n", omp_get_max_threads());
#endif

	Progress p(max, display_progress); // create the progress monitor
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
	for (int i = 0; i < max; ++i) {
		if ( p.increment() ) { // the only way to exit an OpenMP loop
			your_long_computation(nb);
		}
	}
}

RcppExport SEXP test_sequential_wrapper2(SEXP __max, SEXP __nb, SEXP __display_progress) {
	test_sequential2(as<unsigned long>(__max), as<int>(__nb), as<bool>(__display_progress));
	return R_NilValue;
}

RcppExport SEXP test_multithreaded_wrapper2(SEXP __max, SEXP __nb, SEXP __threads, SEXP __display_progress) {
	test_multithreaded_omp2(as<unsigned long>(__max), as<int>(__nb), as<int>(__threads), as<bool>(__display_progress));
	return R_NilValue;
}
