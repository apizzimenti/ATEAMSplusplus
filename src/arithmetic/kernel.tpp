
#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/arithmetic/kernel.h"

#include <SparseRREF/sparse_mat.h>

using namespace std;
using namespace ATEAMS;


template <typename T>
SparseVector<T> arithmetic::randomLinearCombination(
	SparseMatrix<T> kernel,
	const Field& F,
	std::uniform_int_distribution<int>& intuniform,
	std::mt19937& RNG
) {
	// Create a vector of random coefficients, multiply the kernel on the right, done.
	SparseVector<T> rcoeffs;
	for (int i=0; i < kernel.ncol; i++) rcoeffs.push_back((index_t)i, (T)intuniform(RNG));
	rcoeffs.compress();

	return SparseRREF::sparse_mat_dot_sparse_vec<T,index_t>(kernel, rcoeffs, F);
}


template <typename T>
SparseVector<T> arithmetic::submatrixKernelSample(
	SparseMatrix<T> coboundary,
	const Field& F,
	set<size_t> exclude,
	std::uniform_int_distribution<int>& intuniform,
	std::mt19937& RNG,
	arithmetic::ThreadOptions& options,
	bool DEBUG
) {
	// Shrink the coboundary matrix to the appropriate size.
	for (auto i : exclude) coboundary[i].zero();
	coboundary.clear_zero_row(); // TODO check the performance of this!
	coboundary.compress();

	vector<vector<SparsePivot>> pivots;

	// RREF (i.e. find the pivots for) the subcoboundary matrix.
	if (DEBUG) cerr << std::format("RREFing {}x{} coboundary matrix", coboundary.nrow, coboundary.ncol) << endl;
	pivots = SparseRREF::sparse_mat_rref<T,index_t>(coboundary, F, options.opt);

	// Find a kernel for the subcoboundary matrix.
	if (DEBUG) cerr << std::format("computing kernel of {}x{} coboundary matrix", coboundary.nrow, coboundary.ncol) << endl;
	SparseMatrix<T> kernel = SparseRREF::sparse_mat_rref_kernel<T,index_t>(coboundary, pivots, F, options.opt);

	// Return a random linear combination of the columns.
	return arithmetic::randomLinearCombination<T>(kernel, F, intuniform, RNG);
}
