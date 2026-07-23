
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_H
#define ATEAMS_TOPOLOGY_PERSISTENCE_H

#include "ATEAMS++/complexes/Complex.h"

#include <random>
#include <vector>
#include <set>

namespace ATEAMS::topology {
	/**
	 * @brief Computes the persistent homology of a complex using \f$\Z/2\Z\f$
	 * coefficients. See the implementation of @ref ATEAMS::models::Bernoulli.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex (Pointer to) a @ref ATEAMS::complexes::Complex with a flat boundary matrix.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix. For example, if we are doing
	 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
	 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
	 * @param dimension The percolation dimension.
	 * 
	 * @returns A vector of percolation times.
	 */
	template <typename RingLike>
	std::vector<int> PHATPersistence(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		int dimension
	);

	/**
	 * @brief Computes the persistent homology of a complex over \f$\Z/p\Z\f$
	 * coefficients, where \f$p\f$ is prime. Implements the "twist" algorithm
	 * from Chen and Kerber. See the implementation of @ref ATEAMS::models::InvadedCluster
	 * or @ref ATEAMS::models::Invasion.
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
	std::vector<int> twistPersistence(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	);

	/**
	 * @brief Computes the persistent homology of a complex over \f$\Z/p\Z\f$
	 * coefficients, where \f$p\f$ is prime. Implements the standard persistence
	 * algorithm of Edelsbrunner, Harer, and Zomorodian.	
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
	std::vector<int> standardPersistence(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	);

	/**
	 * @brief Dispatch method for computing persistence; calls @ref ATEAMS::topology::PHATPersistence
	 * or @ref ATEAMS::topology::twistPersistence depending on the characteristic
	 * of @p R.	
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix. For example, if we are doing
	 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
	 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
	 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @param dimension The percolation dimension.
	 * 
	 * @returns A vector of percolation times.
	 */
	template <typename RingLike>
	std::vector<int> persistence(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	);

	/** @cond */
	// Expose this for testing, but not for general use.
	template <typename RingLike>
	SparseMatrix<RingLike> reindexSparseBoundaryMatrix(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		int dimension
	);
	/** @endcond */
}

#include "ATEAMS++/topology/persistence.tpp"

#endif

