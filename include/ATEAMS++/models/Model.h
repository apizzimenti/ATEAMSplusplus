
#ifndef ATEAMS_MODELS_MODEL_H
#define ATEAMS_MODELS_MODEL_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/complexes/Complex.h"
#include "ATEAMS++/arithmetic/options.h"

namespace ATEAMS {
	namespace models {
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
			int dimension;
			bool DEBUG = false;
		};

		/**
		 * @struct ModelState
		 * @brief Maintains state for this Model.
		 */
		struct ModelState {};

		/**
		 * @class Model
		 * @brief Abstract (template) class for various Models.
		 * 
		 * @var Model::complex
		 * Pointer to a @ref ATEAMS::Complex.
		 * 
		 * @var Model::parameters
		 * Parameters for the model.
		 * 
		 * @var Model::state
		 * Manages state for the model.
		 */
		class Model {
			public:
				ATEAMS::complexes::Complex* complex;
				ModelParameters parameters;
				ModelState state;

				/**
				 * @brief Draw a sample from the model.
				 * 
				 * @param t When simulated by a Markov chain, \f$t\f$ is the current
				 * 	time-step.
				 * @param options Options for the multithreaded computing environment.
				 * 
				 * @returns A sample from the model.
				 */
				virtual ZpVector sample(int t, ATEAMS::arithmetic::ThreadOptions& options) = 0;

				/** @brief Initializes the state as determined by the model. */
				virtual void initialize() = 0;

				/**
				 * @brief Initializes the state as determined by the user.
				 * 
				 * @param c Initial state.
				 */
				virtual void initialize(ZpVector c) = 0;
		};
	}
}

#endif
