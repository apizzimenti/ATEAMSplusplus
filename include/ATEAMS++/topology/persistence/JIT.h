
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_JIT_H
#define ATEAMS_TOPOLOGY_PERSISTENCE_JIT_H

#include "ATEAMS++/complexes/Complex.h"
#include "ATEAMS++/topology/persistence/reduction.h"

#include <vector>
#include <set>

namespace ATEAMS::topology::persistence {
	/**
	 * @brief Computes the persistent homology of a complex over \f$\Z/p\Z\f$
	 * coefficients, where \f$p\f$ is prime. Implements a parallelized version
	 * of the standard persistence algorithm in line with that of Edelsbrunner
	 * and Harer, incorporating a just-in-time (JIT) version of the twist
	 * optimization. Base implementation to which all overloads forward.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex with a flat boundary matrix.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix. For example, if we are doing
	 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
	 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
	 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @param options Parallel computing environment options.
	 * @param reindexingPolicy How to reindex the full boundary matrix.
	 * @param traversalPolicy How to traverse matrix blocks.
	 * @param reportingPolicy How to report essential cycles.
	 * 
	 * @returns A vector of percolation times, according to @p reportingPolicy.
	 */
	template <typename RingLike>
	std::vector<int> JIT(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options,
		policies::ReindexingPolicy<RingLike> reindexingPolicy,
		policies::TraversalPolicy<RingLike> traversalPolicy,
		policies::ReportingPolicy<RingLike> reportingPolicy
	);

	/**
	 * @brief Computes the persistent homology of a complex over \f$\Z/p\Z\f$
	 * coefficients, where \f$p\f$ is prime. Implements a parallelized version
	 * of the standard persistence algorithm in line with that of Edelsbrunner
	 * and Harer, incorporating a just-in-time (JIT) version of the twist
	 * optimization. Convenience overload for computing persistence of
	 * a single dimension.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex with a flat boundary matrix.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix. For example, if we are doing
	 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
	 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
	 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @param dimension The percolation dimension.
	 * @param options Parallel computing environment options.
	 * 
	 * @returns A vector of percolation times.
	 */
	template <typename RingLike>
	std::vector<int> JIT(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	);


	/**
	 * @brief Computes the persistent homology of a complex over \f$\Z/p\Z\f$
	 * coefficients, where \f$p\f$ is prime. Implements a parallelized version
	 * of the standard persistence algorithm in line with that of Edelsbrunner
	 * and Harer, incorporating a just-in-time (JIT) version of the twist
	 * optimization. Convenience overload.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex with a flat boundary matrix.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix. For example, if we are doing
	 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
	 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
	 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @param options Parallel computing environment options.
	 * 
	 * @returns A vector of percolation times.
	 */
	template <typename RingLike>
	std::vector<int> JIT(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options
	);
}

#include "ATEAMS++/topology/persistence/JIT.tpp"

#endif

