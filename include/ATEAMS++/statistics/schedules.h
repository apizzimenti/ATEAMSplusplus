
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
	}
}


#endif


