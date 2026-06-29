
#ifndef ATEAMS_STATISTICS_TEMPERATURE_H
#define ATEAMS_STATISTICS_TEMPERATURE_H

#include <cmath>

namespace ATEAMS {
	namespace statistics {
		/**
		 * @brief Determines the self-dual probability for the PLGT/PRCM.
		 * 
		 * @param p Order of the finite field for PLGT.
		 * 
		 * @returns A lambda function that always returns the self-dual point
		 * 	\f$\ln(1-\nicefrac{\sqrt p}{1+\sqrt p})\f$.
		 */
		inline std::function<double(int)> selfdual(int p) {
			double c = std::log(
				1-(std::sqrt(p)/(1+std::sqrt(p)))
			);

			return [c](int t) { return c; };
		}

		/**
		 * @brief Constant-temperature.
		 * @param C Temperature to return. If you want to see whether this temperature
		 * 	is valid (i.e. \f$0 \leq 1-e^{-C} \leq 1\f$), please set the `DEBUG`
		 * 	model parameter (e.g. in @ref ATEAMS::models::SwendsenWangParameters)
		 *  to `true`. If your provided temperature doesn't work, an assertion will
		 * 	fail.
		 * @returns Lambda function that always returns \f$C\f$.
		 */
		inline std::function<double(int)> constant(double C) {
			return [C](int t) { return C; };
		}

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


