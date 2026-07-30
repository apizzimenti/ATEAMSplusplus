
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_H
#define ATEAMS_TOPOLOGY_PERSISTENCE_H

#include "ATEAMS++/complexes/Complex.h"

#include "ATEAMS++/topology/persistence/JIT.h"
#include "ATEAMS++/topology/persistence/parallel.h"
#include "ATEAMS++/topology/persistence/PHAT.h"
#include "ATEAMS++/topology/persistence/split.h"
#include "ATEAMS++/topology/persistence/stagger.h"
#include "ATEAMS++/topology/persistence/standard.h"
#include "ATEAMS++/topology/persistence/twist.h"

#include <vector>
#include <set>


/** @brief Persistent homology routines. */
namespace ATEAMS::topology::persistence {
	/**
	 * @brief Dispatch method for computing persistence; dispatches to 
	 * @ref ATEAMS::topology::persistence::PHAT, @ref ATEAMS::topology::persistence::twist,
	 * or @ref ATEAMS::topology::persistence::JIT depending on the characteristic
	 * of @p R and parallelization options.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
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
	inline std::vector<int> persistence(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		std::vector<int> essential;

		if (R->characteristic == 2) essential = PHAT<RingLike>(complex, filtration);
		else {
			if (options.parallel->enabled) essential = JIT<RingLike>(complex, filtration, R, options);
			else essential = twist<RingLike>(complex, filtration, R, options);
		}

		std::sort(essential.begin(), essential.end());
		return essential;
	};

	/**
	 * @brief Dispatch method for computing persistence; dispatches to 
	 * @ref ATEAMS::topology::persistence::PHAT, @ref ATEAMS::topology::persistence::twist,
	 * or @ref ATEAMS::topology::persistence::JIT depending on the characteristic
	 * of @p R and parallelization options. Convenience overload for computing
	 * the persistence of a single dimension.
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
	 * @param options Options for the multithreaded computing environment.
	 * 
	 * @returns A vector of percolation times.
	 */
	template <typename RingLike>
	inline std::vector<int> persistence(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		std::vector<int> essential;

		if (R->characteristic == 2) essential = PHAT<RingLike>(complex, filtration, dimension);
		else {
			if (options.parallel->enabled) essential = JIT<RingLike>(complex, filtration, R, dimension, options);
			else essential = twist<RingLike>(complex, filtration, R, dimension, options);
		}

		std::sort(essential.begin(), essential.end());
		return essential;
	};
}

#endif

