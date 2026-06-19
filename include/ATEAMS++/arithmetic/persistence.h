
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
		std::vector<int> PHATPersistence(
			ATEAMS::complexes::Complex* complex,
			std::vector<int> filtration,
			int dimension
		);

		/**
		 * @brief Computes the persistent homology of a complex over \f$\Z/q\Z\f$
		 * 	coefficients, where \f$q\f$ is prime.
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
		std::vector<int> TwistPersistence(
			ATEAMS::complexes::Complex* complex,
			std::vector<int> filtration,
			Zp F,
			int dimension
		);

		/**
		 * @brief A convenience function for deciding when to stop invading
		 * 	(i.e. when to re-sample the spins over the included edges).
		 * 
		 * @param ranks A vector of integers with values in \f$[1,\rank(H_d(X))]\f$,
		 * 	where \f$\rank(H_d(X))\f$ is the rank of the \f$d\f$th homology group
		 * 	of \f$X\f$.
		 * 
		 * @returns A function that draws a value from `ranks` at uniform random,
		 * 	then subtracts \f$1\f$.
		 */
		inline std::function<int(int)> stopInvadingAt(std::vector<int> ranks) {
			return [ranks](int t) {
				std::random_device rd;
				std::mt19937 RNG(rd());
				std::uniform_int_distribution<int> uniform(0,ranks.size()-1);

				return ranks[uniform(RNG)]-1;
			};
		}
	}
}

#endif

