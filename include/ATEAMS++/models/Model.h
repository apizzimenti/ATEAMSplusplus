
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
	 * @var ModelParameters::dimension
	 * Dimension of the cells out of which samples are built.
	 * 
	 * @var ModelParameters::DEBUG
	 * Are we debugging the Model?
	 */
	struct ModelParameters {
		// Used by: SwendsenWang, InvadedCluster, Invasion, Glauber
		int field;

		// Used by: SwendsenWang, Glauber
		std::function<double(int)> temperatureFunction;

		// Universal.
		int dimension;
		bool DEBUG = false;
	};

	/**
	 * @struct ModelState
	 * @brief Maintains state for this Model.
	 */
	template <typename Container>
	struct ModelState {
		// Used by: SwendsenWang, InvadedCluster,
		Container cochain;

		// Used by: SwendsenWang, InvadedCluster, Invasion, Bernoulli
		std::vector<int> includes;

		// Used by: Bernoulli, Invasion
		std::vector<int> essential;

		// Used by: Bernoulli
		int rank;

		// Used by: Glauber
		int energy;
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
	 * @var Model::state
	 * Manages state for the model.
	 * 
	 * @var Model::kind
	 * Model name.
	 */
	template <typename T, template <typename> typename ContainerType>
	class Model {
		public:
			ATEAMS::complexes::Complex<T>* complex;
			ModelParameters parameters;
			ModelState<ContainerType<T>> state;
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
			virtual ModelState<ContainerType<T>> sample(int t, ATEAMS::arithmetic::ThreadOptions& options) = 0;

			/** @brief Initializes the state as determined by the model. */
			virtual void initialize() = 0;

			/**
			 * @brief Initializes the state as determined by the user.
			 * 
			 * @param c Initial state.
			 */
			virtual void initialize(ContainerType<T> c) = 0;
	};
}

#endif
