
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
		 * @tparam RingLike A coefficient @ref Ring type, like @ref Zp or @ref Q.
		 * @tparam VectorLike Storage unit compatible with @p RingLike, like
	 	 *  @ref SparseVector or @ref DenseVector.
		 * 
		 * **Example Usage**
		 * 
		 * @code
		 * using namespace ATEAMS;
		 * 
		 * using Model = models::Bernoulli;
		 * using Parameters = models::ModelParameters;

		 * using Complex = complexes::Cubical<Model::RingType>;
		 * using State = models::ModelState<Model::RingType,Model::VectorType>;
		 * using Chain = statistics::Chain<Model::RingType,Model::VectorType>;
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
		 * for (State state : M.simulate()) {
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
		 * 	@brief Multithreaded computing options @ref arithmetic::ThreadOptions.
		 * 
		 * @var Chain::steps
		 * 	@brief Number of iterations.
		 */
		template <typename ModelType>
		class Chain {
			public:
				// models::Model<typename ModelType::RingType,typename ModelType::VectorType>* model;
				ModelType* model;
				typename ModelType::State state;
				arithmetic::ThreadOptions options;

				int steps;

				/**
				 * @brief Constructor; uses default @ref arithmetic::ThreadOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 */
				Chain(ModelType* model, int steps) {
					this->model = model;
					this->steps = steps;

					arithmetic::ThreadOptions options;
					this->options = options;

					typename ModelType::State state;
					this->state = state;
				};

				/**
				 * @brief Constructor; uses user-provided @ref ATEAMS::arithmetic::ThreadOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 * @param options User-provided compute options.
				 */
				Chain(ModelType* model, int steps, arithmetic::ThreadOptions options) {
					this->model = model;
					this->steps = steps;
					this->options = options;

					typename ModelType::State state;
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
				std::generator<typename ModelType::State> simulate() {
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

