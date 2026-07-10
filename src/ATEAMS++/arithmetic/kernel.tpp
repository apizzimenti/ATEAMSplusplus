
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

namespace ATEAMS::arithmetic {
	template <typename RingLike>
	SparseVector<RingLike> randomLinearCombination(
		SparseMatrix<RingLike> kernel,
		Ring* R,
		std::uniform_int_distribution<int>& intuniform,
		std::mt19937& RNG
	) {
		// Create a vector of random coefficients, multiply the kernel on the right, done.
		SparseVector<RingLike> rcoeffs;
		for (int i=0; i < kernel.ncol; i++) rcoeffs.push_back((INDEX)i, (typename RingLike::dtype)intuniform(RNG));
		rcoeffs.compress();

		return SparseRightMultiplication<RingLike>(kernel, rcoeffs, R);
	}


	template <typename RingLike>
	SparseVector<RingLike> submatrixKernelSample(
		SparseMatrix<RingLike> coboundary,
		Ring* R,
		set<size_t> exclude,
		std::uniform_int_distribution<int>& intuniform,
		std::mt19937& RNG,
		ThreadOptions& options,
		bool DEBUG
	) {
		// Shrink the coboundary matrix to the appropriate size.
		for (auto i : exclude) coboundary[i].zero();
		coboundary.clear_zero_row(); // TODO check the performance of this!
		coboundary.compress();

		SparseMatrix<RingLike> kernel;

		if (DEBUG) {
			// RREF (i.e. find the pivots for) the subcoboundary matrix.
			cerr << std::format("RREFing {}x{} coboundary matrix", coboundary.nrow, coboundary.ncol) << endl;
			SparsePivots pivots = SparseMatrixRREF<RingLike>(coboundary, R, options);

			// Find a kernel for the subcoboundary matrix.
			cerr << std::format("computing kernel of {}x{} coboundary matrix", coboundary.nrow, coboundary.ncol) << endl;
			kernel = SparseMatrixRREFKernel<RingLike>(coboundary, R, pivots, options);
		} else {
			kernel = SparseMatrixKernel<RingLike>(coboundary, R, options);
		}

		// Return a random linear combination of the columns.
		return randomLinearCombination<RingLike>(kernel, R, intuniform, RNG);
	}
}

#endif
