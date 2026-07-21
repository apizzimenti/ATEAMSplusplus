
#ifndef ATEAMS_STATISTICS_CHAIN_H
#define ATEAMS_STATISTICS_CHAIN_H

#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/arithmetic.h"

#include <generator>


namespace ATEAMS {
	namespace statistics {
		/**
		 * @class Chain
		 * @brief Simulates the Markov chains specified by each @ref models::Model.
		 * 
		 * @tparam ModelType A @ref ATEAMS::models::Model type, like @ref ATEAMS::models::Bernoulli.
		 * 
		 * **Example Usage**
		 * 
		 * The below code draws 1,000 samples from plaquette Bernoulli percolation
		 * on \f$\T^4_6\f$, the \f$4\f$-fold torus of scale \f$6\f$.
		 * @code
		 * using namespace ATEAMS;
		 * 
		 * using Model = models::Bernoulli;
		 * using Complex = complexes::Cubical<Model::RingType>;
		 * 
		 * using Parameters = models::ModelParameters;
		 * using Chain = statistics::Chain<Model>;
		 * 
		 * Complex plex({4,4,4,4});
		 * 
		 * Parameters params;
		 * params.p = 0.5;
		 * params.dimension = 2;
		 * 
		 * Model percolation(&plex, parameters);
		 * Chain M(&percolation, 1000);
		 * 
		 * for (Chain::State state : M.simulate()) {
		 * 		<...>
		 * }
		 * @endcode
		 * 
		 * @var Chain::model
		 * 	@brief (Pointer to) a Model.
		 * 
		 * @var Chain::state
		 * 	@brief A @ref models::ModelState.
		 * 
		 * @var Chain::options
		 * 	@brief Multithreaded computing options @ref arithmetic::ComputeOptions.
		 * 
		 * @var Chain::steps
		 * 	@brief Number of iterations.
		 */
		template <typename ModelType>
		class Chain {
			public:
				ModelType* model;
				arithmetic::ComputeOptions<typename ModelType::RingType> options;

				int steps;

				/**
				 * @brief Exposed typename for the model state.
				 */
				using State = typename ModelType::State;
				State state;

				/**
				 * @brief Constructor; uses default @ref arithmetic::ComputeOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 */
				Chain(ModelType* model, int steps) {
					this->model = model;
					this->steps = steps;

					arithmetic::ComputeOptions<typename ModelType::RingType> options;
					this->options = options;

					State state;
					this->state = state;
				};

				/**
				 * @brief Constructor; uses user-provided @ref ATEAMS::arithmetic::ComputeOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 * @param options User-provided compute options.
				 */
				Chain(
					ModelType* model,
					int steps,
					arithmetic::ComputeOptions<typename ModelType::RingType> options
				) {
					this->model = model;
					this->steps = steps;
					this->options = options;

					State state;
					this->state = state;
				};

				/**
				 * @brief Model-generic iterator.
				 * 
				 * @code
				 * using Model = models::Bernoulli;
				 * using State = models::ModelState<Model::RingType,Model::VectorType>;
				 * 
				 * for (State state : M.simulate()) {
				 *     <do whatever>
				 * }
				 * @endcode
				 * 
				 * @returns A `std::generator`.
				 */
				std::generator<State> simulate() {
					std::thread listener = options.spinUp();
					this->model->initialize(this->state);

					for (int t=0; t < this->steps; t++) {
						this->model->sample(t, this->state, this->options);
						co_yield this->state;
					}

					// Spin down multithreading environment.
					options.spinDown(&listener);
				}
		};
	}
}


#endif

