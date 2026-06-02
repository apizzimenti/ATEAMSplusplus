
#ifndef ATEAMS_STATISTICS_TEMPERATURE_H
#define ATEAMS_STATISTICS_TEMPERATURE_H

#include <cmath>

namespace ATEAMS {

	/**
	 * Implements a critical-temperature function for the PLGT.
	 */
	inline std::function<double(int)> critical(int p) {
		double c = std::log(
			1-(std::sqrt(p)/(1+std::sqrt(p)))
		);

		return [c](int t) { return c; };
	}

	/**
	 * Implements a constant-temperature function for the PLGT.
	 */
	inline std::function<double(int)> constant(double C) {
		return [C](int t) { return C; };
	}
}


#endif


