
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_PHAT_H
#define ATEAMS_TOPOLOGY_PERSISTENCE_PHAT_H

#include "ATEAMS++/complexes/Complex.h"
#include "ATEAMS++/topology/persistence/reduction.h"

#include <random>
#include <vector>
#include <set>


/** @brief Persistent homology routines. */
namespace ATEAMS::topology::persistence {
	/**
	 * @brief Computes the persistent homology of the @p dimension-skeleton of
	 * @p complex with \f$\Z/2\Z\f$ coefficients using <a href="https://tinyurl.com/223wxnye">PHAT</a>.
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
	std::vector<int> PHAT(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		int dimension
	);

	/**
	 * @brief Computes the persistent homology of a complex using \f$\Z/2\Z\f$
	 * coefficients using <a href="https://tinyurl.com/223wxnye">PHAT</a>. See
	 * the implementation of @ref ATEAMS::models::Bernoulli.
	 * 
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex (Pointer to) a @ref ATEAMS::complexes::Complex with a flat boundary matrix.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix.
	 * 
	 * @returns A vector of percolation times.
	 */
	template <typename RingLike>
	std::vector<int> PHAT(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration
	);
}

#include "ATEAMS++/topology/persistence/PHAT.tpp"

#endif

