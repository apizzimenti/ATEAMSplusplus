
#ifndef ATEAMS_ARITHMETIC_KERNEL_H
#define ATEAMS_ARITHMETIC_KERNEL_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/arithmetic/options.h"

#include <vector>
#include <random>
#include <set>

namespace ATEAMS::arithmetic {
	/**
	 * @brief A uniform random linear combination of the rows of \f$K\f$,
	 * typically the kernel of a reduced coboundary matrix \f$\delta^{d-1}\f$.
	 * See the implementation of @ref ATEAMS::arithmetic::submatrixKernelSample.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param K Kernel.
	 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @param intuniform Uniform distribution over integers \f$[0, p)\f$.
	 * @param RNG Mersenne twister PRNG.
	 * 
	 * @return (Sparse) uniform vector in \f$\ker(\delta^{d-1})\f$.
	 */
	template <typename RingLike>
	SparseVector<RingLike> randomLinearCombination(
		SparseMatrix<RingLike> K,
		Ring* R,
		std::uniform_int_distribution<int>& intuniform,
		std::mt19937& RNG
	);
	
	/**
	 * @brief Draws a cochain \f$f\f$ at uniform random from \f$\delta^{d-1}\f$
	 * such that \f$f\f$ is a cocycle on the \f$d\f$-cells corresponding to
	 * the rows not in `exclude`. See the implementations of
	 * @ref ATEAMS::models::SwendsenWang or @ref ATEAMS::models::InvadedCluster
	 * for usage.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param coboundary Coboundary matrix \f$\delta^{d-1}\f$. See @ref ATEAMS::complexes::BoundaryData.
	 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @param exclude Row indices to exclude.
	 * @param intuniform Uniform distribution over integers \f$[0, p)\f$.
	 * @param RNG Mersenne twister PRNG.
	 * @param options Reference to a @ref arithmetic::ComputeOptions instance.
	 * @param DEBUG _(Optional, default `false`)_ Are we debugging the arithmetic?
	 * 
	 * @return (Sparse) uniform vector in \f$\ker(\delta^{d-1})\f$.
	 */
	template <typename RingLike>
	SparseVector<RingLike> submatrixKernelSample(
		SparseMatrix<RingLike> coboundary,
		Ring* R,
		std::set<size_t> exclude,
		std::uniform_int_distribution<int>& intuniform,
		std::mt19937& RNG,
		ComputeOptions<RingLike>& options,
		bool DEBUG=false
	);
}

#include "ATEAMS++/arithmetic/kernel.tpp"

#endif
