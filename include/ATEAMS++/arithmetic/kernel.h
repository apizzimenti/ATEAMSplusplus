
#ifndef ATEAMS_ARITHMETIC_KERNEL_H
#define ATEAMS_ARITHMETIC_KERNEL_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/arithmetic/options.h"

#include <vector>
#include <random>
#include <set>

namespace ATEAMS {
	namespace arithmetic {
		/**
		 * @brief Generates a uniform random linear combination of the basis for
		 * the (reduced) coboundary matrix \f$\delta^{d-1}\f$.
		 * 
		 * @param kernel Kernel of the reduced coboundary matrix.
		 * @param F Field \f$ \Z / p \Z \f$.
		 * @param intuniform Uniform distribution over integers \f$[0, p)\f$.
		 * @param RNG Mersenne twister PRNG.
		 * 
		 * @return (Sparse) uniform vector in \f$\ker(\delta^{d-1})\f$.
		 */
		template <typename T>
		SparseVector<T> randomLinearCombination(
			SparseMatrix<T> kernel,
			const Field& F,
			std::uniform_int_distribution<int>& intuniform,
			std::mt19937& RNG
		);
		
		/**
		 * @brief Draws a cochain \f$f\f$ at uniform random from \f$\delta^{d-1}\f$
		 * such that \f$f\f$ is a cocycle on the \f$d\f$-cells corresponding to
		 * the rows not in `exclude`.
		 * 
		 * @param coboundary Coboundary matrix \f$\delta^{d-1}\f$.
		 * @param F Field \f$ \Z / p \Z \f$.
		 * @param exclude Row indices to exclude.
		 * @param intuniform Uniform distribution over integers \f$[0, p)\f$.
		 * @param RNG Mersenne twister PRNG.
		 * @param options Reference to a @ref ThreadOptions.
		 * 	instance.
		 * @param DEBUG _(Optional)_ Are we debugging the arithmetic?
		 * 
		 * @return (Sparse) uniform vector in \f$\ker(\delta^{d-1})\f$.
		 */
		template <typename T>
		SparseVector<T> submatrixKernelSample(
			SparseMatrix<T> coboundary,
			const Field& F,
			std::set<size_t> exclude,
			std::uniform_int_distribution<int>& intuniform,
			std::mt19937& RNG,
			ATEAMS::arithmetic::ThreadOptions& options,
			bool DEBUG=false
		);
	}
}

#include "ATEAMS++/arithmetic/kernel.tpp"

#endif
