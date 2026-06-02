
#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/arithmetic/kernel.h"
#include "libraries/SparseRREF/sparse_mat.h"

using namespace std;

/**
 * Creates a random linear combination of the columns of the kernel passed to it.
 */
ZpVector ATEAMS::randomLinearCombination(
	ZpMatrix kernel,
	const Zp& F,
	std::uniform_int_distribution<int>& intuniform,
	std::mt19937& RNG
) {
	// Create a vector of random coefficients, multiply the kernel on the right, done.
	ZpVector rcoeffs;
	for (int i=0; i < kernel.ncol; i++) rcoeffs.push_back((index_t)i, (data_t)intuniform(RNG));
	rcoeffs.compress();

	return SparseRREF::sparse_mat_dot_sparse_vec<data_t, index_t>(kernel, rcoeffs, F);
}


/**
 * Given a coboundary matrix and a list of rows to exclude, samples a uniform
 * random element of the kernel of the coboundary matrix with those rows excluded.
 */
ZpVector ATEAMS::submatrixKernelSample(
	ZpMatrix coboundary,
	const Zp& F,
	set<size_t> exclude,
	std::uniform_int_distribution<int>& intuniform,
	std::mt19937& RNG,
	bool _DEBUG
) {
	// FLINT memory stuff.
	Flint::set_memory_functions();
	SparseRREF::rref_option_t opt;
	opt->pool.reset();
	opt->method = 0;

	// Listen for the "abort" event.
	std::thread listener(key_listener, ref(opt->abort));

	// Shrink the coboundary matrix to the appropriate size.
	for (auto i : exclude) coboundary[i].zero();
	coboundary.clear_zero_row(); // TODO check the performance of this!
	coboundary.compress();

	vector<vector<ZpPivot>> pivots;

	// RREF (i.e. find the pivots for) the subcoboundary matrix.
	if (_DEBUG) cerr << std::format("RREFing {}x{} coboundary matrix", coboundary.nrow, coboundary.ncol) << endl;
	pivots = SparseRREF::sparse_mat_rref<data_t, index_t>(coboundary, F, opt);

	// Find a kernel for the subcoboundary matrix.
	if (_DEBUG) cerr << std::format("computing kernel of {}x{} coboundary matrix", coboundary.nrow, coboundary.ncol) << endl;
	ZpMatrix kernel = SparseRREF::sparse_mat_rref_kernel<data_t, index_t>(coboundary, pivots, F, opt);

	// Spin down FLINT memory stuff.
	opt->abort = true;
	Flint::clear_cache();
	listener.join();

	// Return a random linear combination of the columns.
	return ATEAMS::randomLinearCombination(kernel, F, intuniform, RNG);
}
