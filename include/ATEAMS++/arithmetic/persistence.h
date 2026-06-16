
#ifndef ATEAMS_ARITHMETIC_PERSISTENCE_H
#define ATEAMS_ARITHMETIC_PERSISTENCE_H

#include "ATEAMS++/complexes/complexes.h"

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
		std::vector<int> PHATPersistence(
			ATEAMS::complexes::Complex* complex,
			std::vector<int> filtration,
			int dimension
		);
	}
}

#endif

