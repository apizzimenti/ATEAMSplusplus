
#ifndef ATEAMS_MODELS_MODELPARAMETERS_H
#define ATEAMS_MODELS_MODELPARAMETERS_H

#include "ATEAMS++/common.h"

namespace ATEAMS::models {
	/**
	 * @struct ModelParameters
	 * @brief Typical parameters required for simulating a Model.
	 * 
	 * @var ModelParameters::coefficients
	 * 	@brief Coefficient ring.
	 * 
	 * @var ModelParameters::dimension
	 * 	Percolation subcomplex dimension. Used by every @ref Model.
	 * 
	 * @var ModelParameters::p
	 * 	Bernoulli trial density. Used by @ref Bernoulli.
	 * 
	 * @var ModelParameters::temperatureFunction
	 * 	Function that specifies the (inverse) temperature parameter at the current
	 * 	time-step. Used by @ref Glauber, @ref SwendsenWang.
	 * 
	 * @var ModelParameters::stoppingFunction
	 * 	Function that specifies the number of giant cycles encountered before
	 * 	re-sampling spins. Used by @ref InvadedCluster, @ref Invasion.
	 * 
	 * @var ModelParameters::DEBUG
	 * Are we debugging the @ref Model?
	 */
	struct ModelParameters {
		// Used by: SwendsenWang, InvadedCluster, Invasion, Glauber
		Ring* coefficients;

		// Used by: Bernoulli
		float p;

		// Used by: SwendsenWang, Glauber
		std::function<double(int)> temperatureFunction;

		// Used by: InvadedCluster, Invasion
		std::function<int(int)> stoppingFunction;

		// Universal.
		int dimension;
		bool DEBUG = false;
	};
}

#endif
