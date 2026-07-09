
#ifndef ATEAMS_ARITHMETIC_KERNEL_T
#define ATEAMS_ARITHMETIC_KERNEL_T

#ifndef ATEAMS_ARITHMETIC_KERNEL_H
#error __FILE__ should only be included from arithmetic/kernel.h.
#endif

#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/arithmetic/kernel.h"

#include <SparseRREF/sparse_mat.h>

using namespace std;
using namespace ATEAMS;


template <typename RingLike>
SparseVector<RingLike> arithmetic::randomLinearCombination(
	SparseMatrix<RingLike> kernel,
	Ring* R,
	std::uniform_int_distribution<int>& intuniform,
	std::mt19937& RNG
) {
	// Create a vector of random coefficients, multiply the kernel on the right, done.
	SparseVector<RingLike> rcoeffs;
	for (int i=0; i < kernel.ncol; i++) rcoeffs.push_back((INDEX)i, (typename RingLike::dtype)intuniform(RNG));
	rcoeffs.compress();

	return SparseRREF::sparse_mat_dot_sparse_vec<typename RingLike::dtype,INDEX>(kernel, rcoeffs, R->ring);
}


template <typename RingLike>
SparseVector<RingLike> arithmetic::submatrixKernelSample(
	SparseMatrix<RingLike> coboundary,
	Ring* R,
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
	pivots = SparseRREF::sparse_mat_rref<typename RingLike::dtype,INDEX>(coboundary, R->ring, options.opt);

	// Find a kernel for the subcoboundary matrix.
	if (DEBUG) cerr << std::format("computing kernel of {}x{} coboundary matrix", coboundary.nrow, coboundary.ncol) << endl;
	SparseMatrix<RingLike> kernel = SparseRREF::sparse_mat_rref_kernel<typename RingLike::dtype,INDEX>(coboundary, pivots, R->ring, options.opt);

	// Return a random linear combination of the columns.
	return arithmetic::randomLinearCombination<RingLike>(kernel, R, intuniform, RNG);
}

#endif
