
#ifndef ATEAMS_STATISTICS_CHAIN_H
#define ATEAMS_STATISTICS_CHAIN_H

#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <generator>


namespace ATEAMS {
	namespace statistics {

		/**
		 * @class Chain
		 * @brief Simulates the Markov chains specified by each @ref ATEAMS::models::Model.
		 * 
		 * @var Chain::model
		 * 	(Pointer to) a Model.
		 * 
		 * @var Chain::state
		 * 	(Pointer to) Model state.
		 * 
		 * @var Chain::options
		 * 	Multithreaded computing options @ref ATEAMS::arithmetic::ThreadOptions.
		 * 
		 * @var Chain::steps
		 * 	Number of iterations.
		 */
		template <typename ModelType>
		class Chain {
			public:
				ModelType* model;
				ATEAMS::models::ModelState* state;
				ATEAMS::arithmetic::ThreadOptions options;

				int steps;

				/**
				 * @brief Constructor; uses default @ref ATEAMS::arithmetic::ThreadOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 */
				Chain(ModelType* model, int steps) {
					this->model = model;
					this->steps = steps;

					ATEAMS::arithmetic::ThreadOptions options;
					this->options = options;
				};

				/**
				 * @brief Constructor; use user-provided @ref ATEAMS::arithmetic::ThreadOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 */
				Chain(ModelType* model, int steps, ATEAMS::arithmetic::ThreadOptions options) {
					this->model = model;
					this->steps = steps;
					this->options = options;
				};

				/**
				 * @brief Model-generic iterator.
				 * 
				 * @code
				 * 	for (models::ModelState* _state : M.simulate()) {
				 * 		// Requires a cast to access specific members of the state object.
				 * 		models::InvadedClusterState* state = (models::InvadedClusterState*)_state;
				 * 	}
				 * @endcode
				 * 
				 * @returns A `std::generator`.
				 */
				std::generator<ATEAMS::models::ModelState*> simulate() {
					std::thread listener = options.spinUp();

					for (int t=0; t < 1000; t++) {
						this->model->sample(t, this->options);
						this->state = &this->model->state;
						co_yield this->state;
					}

					// Spin down multithreading environment.
					options.spinDown(&listener);
				}

				/**
				 * @brief Templated iterator.
				 * 
				 * @code
				 * 	for (models::InvadedClusterState* _state : M.run<models::InvadedClusterState*>()) {
				 * 		// Doesn't require a cast.
				 * 		...
				 * 	}
				 * @endcode
				 * 
				 * @returns A `std::generator`.
				 */
				template <typename State>
				std::generator<State*> simulate() {
					std::thread listener = options.spinUp();

					for (int t=0; t < 1000; t++) {
						this->model->sample(t, this->options);
						this->state = &this->model->state;
						co_yield (State*)this->state;
					}

					// Spin down multithreading environment.
					options.spinDown(&listener);
				}
		};
	}
}


#endif

