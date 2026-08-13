
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_PARALLEL_H
#define ATEAMS_TOPOLOGY_PERSISTENCE_PARALLEL_H

#include "ATEAMS++/complexes/Complex.h"
#include "ATEAMS++/topology/persistence/reduction.h"

#include <vector>
#include <set>

namespace ATEAMS::topology::persistence {
	/**
	 * @brief Computes the persistent homology of a complex over \f$\Z/p\Z\f$
	 * coefficients, where \f$p\f$ is prime. Implements a parallelized version
	 * of the standard persistence algorithm in line with that of Edelsbrunner
	 * and Harer. Base implementation to which all overloads forward.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex with a flat boundary matrix.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix. For example, if we are doing
	 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
	 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
	 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @param resources Parallel computing environment resources.
	 * @param reindexingPolicy How to reindex the full boundary matrix.
	 * @param traversalPolicy How to traverse matrix blocks.
	 * @param reportingPolicy How to report essential cycles.
	 * 
	 * @returns A vector of percolation times, according to @p reportingPolicy.
	 */
	template <typename RingLike>
	inline std::vector<int> parallel(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources,
		auto& reindexingPolicy,
		auto& traversalPolicy,
		auto& reportingPolicy
	);


	/**
	 * @brief Computes the persistent homology of a complex over \f$\Z/p\Z\f$
	 * coefficients, where \f$p\f$ is prime. Implements a parallelized version
	 * of the standard persistence algorithm in line with that of Edelsbrunner
	 * and Harer. Overload for computing the persistent homology of a specific
	 * dimension.
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
	 * @param resources Parallel computing environment resources.
	 * 
	 * @returns A vector of percolation times.
	 */
	template <typename RingLike>
	inline std::vector<int> parallel(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeResources<RingLike>& resources
	);

	/**
	 * @brief Computes the persistent homology of a complex over \f$\Z/p\Z\f$
	 * coefficients, where \f$p\f$ is prime. Implements a parallelized version
	 * of the standard persistence algorithm in line with that of Edelsbrunner
	 * and Harer. Convenience overload.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex with a flat boundary matrix.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix. For example, if we are doing
	 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
	 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
	 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @param resources Parallel computing environment resources.
	 * 
	 * @returns A vector of percolation times.
	 */
	template <typename RingLike>
	inline std::vector<int> parallel(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources
	);
}

#include "ATEAMS++/topology/persistence/parallel.tpp"

#endif

