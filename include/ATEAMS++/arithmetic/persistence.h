
#ifndef ATEAMS_ARITHMETIC_PERSISTENCE_H
#define ATEAMS_ARITHMETIC_PERSISTENCE_H

#include "ATEAMS++/complexes/complexes.h"

#include <random>
#include <vector>
#include <set>

namespace ATEAMS {
	namespace arithmetic {
		/**
		 * @brief Computes the persistent homology of a complex using \f$\Z/2\Z\f$
		 * coefficients.
		 * 
		 * @param complex A complex with a flat boundary matrix.
		 * @param filtration A vector that specifies the order in which to add
		 * 	the cells in the flat boundary matrix. For example, if we are doing
		 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
		 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
		 * @param dimension The percolation dimension.
		 * 
		 * @returns A vector of percolation times.
		 */
		template <typename T>
		std::vector<int> PHATPersistence(
			ATEAMS::complexes::Complex<T>* complex,
			std::vector<int> filtration,
			int dimension
		);

		/**
		 * @brief Computes the persistent homology of a complex over \f$\Z/p\Z\f$
		 * 	coefficients, where \f$p\f$ is prime.
		 * 
		 * @param complex A complex with a flat boundary matrix.
		 * @param filtration A vector that specifies the order in which to add
		 * 	the cells in the flat boundary matrix. For example, if we are doing
		 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
		 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
		 * @param F Finite field \f$\Z/p\Z\f$ over which we're doing arithmetic.
		 * @param dimension The percolation dimension.
		 * 
		 * @returns A vector of percolation times.
		 */
		template <typename T>
		std::vector<int> TwistPersistence(
			ATEAMS::complexes::Complex<T>* complex,
			std::vector<int> filtration,
			Field F,
			int dimension
		);

		/** @cond */
		// Expose this for testing, but not for general use.
		template <typename T>
		ATEAMS::SparseMatrix<T> reindexSparseBoundaryMatrix(ATEAMS::complexes::Complex<T>* complex, std::vector<int> filtration, int dimension);
		/** @endcond */

		#include "arithmetic/persistence.tpp"
	}
}

#endif

