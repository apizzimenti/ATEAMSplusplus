
#ifndef ATEAMS_MODELS_MODEL_H
#define ATEAMS_MODELS_MODEL_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/complexes/Complex.h"
#include "ATEAMS++/arithmetic/options.h"

#include <string>

namespace ATEAMS::models {
	/**
	 * @struct ModelParameters
	 * @brief Typical parameters required for simulating a Model.
	 * 
	 * @var ModelParameters::field
	 * 	Field over which computations are performed. Used by every @ref Model.
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
		int field;

		// Used by: Bernoulli
		int p;

		// Used by: SwendsenWang, Glauber
		std::function<double(int)> temperatureFunction;

		// Used by: InvadedCluster, Invasion
		std::function<int(int)> stoppingFunction;

		// Universal.
		int dimension;
		bool DEBUG = false;
	};

	/**
	 * @struct ModelState
	 * @brief Maintains state for this Model.
	 * 
	 * @var ModelState::cochain
	 * 	Current cochain \f$ f_t \f$ in `ContainerType<T>`. Used by
	 * 	@ref InvadedCluster, @ref SwendsenWang.
	 * 
	 * @var ModelState::includes
	 * 	Percolation subcomplex \f$ P_t \f$. Used by @ref Bernoulli, @ref InvadedCluster,
	 * 	@ref Invasion, @ref SwendsenWang.
	 * 
	 * @var ModelState::essential
	 * 	Times at which essential cycles of \f$ P_t \f$ were born. Used by
	 * 	@ref Bernoulli, @ref Invasion.
	 * 
	 * @var ModelState::rank
	 * 	Rank of the \f$d\f$th persistent homology group \f$PH_d(P_t)\f$. Used by
	 * 	@ref Bernoulli.
	 * 
	 * @var ModelState::energy
	 * 	Energy of the current cochain, i.e. \f$\H(f_t)\f$. Used by @ref Glauber.
	 * 
	 * @var ModelState::time
	 * 	Current time-step. Used by every @ref Model.
	 */
	template <typename T, template <typename> typename ContainerType>
	struct ModelState {
		// Used by: SwendsenWang, InvadedCluster,
		ContainerType<T> cochain;

		// Used by: SwendsenWang, InvadedCluster, Invasion, Bernoulli
		std::vector<int> includes;

		// Used by: Bernoulli, Invasion
		std::vector<int> essential;

		// Used by: Bernoulli
		int rank;

		// Used by: Glauber
		int energy;

		// Used by: everyone
		int t;
	};

	/**
	 * @class Model
	 * @brief Abstract (template) class for various Models.
	 * 
	 * @var Model::complex
	 * Pointer to a @ref ATEAMS::complexes::Complex.
	 * 
	 * @var Model::parameters
	 * Parameters for the model.
	 * 
	 * @var Model::kind
	 * Model name.
	 */
	template <typename T, template <typename> typename ContainerType>
	class Model {
		public:
			ATEAMS::complexes::Complex<T>* complex;
			ModelParameters parameters;
			std::string kind;

			/**
			 * @brief Draw a sample from the model.
			 * 
			 * @param t When simulated by a Markov chain, \f$t\f$ is the current
			 * 	time-step.
			 * @param options Options for the multithreaded computing environment.
			 * 
			 * @returns A sample from the model.
			 */
			virtual ModelState<T,ContainerType> sample(int t, ModelState<T,ContainerType>& state, ATEAMS::arithmetic::ThreadOptions& options) = 0;

			/** @brief Initializes the state as determined by the model. */
			virtual ModelState<T,ContainerType> initialize(ModelState<T,ContainerType>& state) = 0;

			/**
			 * @brief Initializes the state as determined by the user.
			 * 
			 * @param c Initial state.
			 */
			virtual ModelState<T,ContainerType> initialize(ContainerType<T> c, ModelState<T,ContainerType>& state) = 0;
	};
}

#endif
