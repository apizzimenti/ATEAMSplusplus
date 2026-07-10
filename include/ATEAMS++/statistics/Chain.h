
#ifndef ATEAMS_STATISTICS_CHAIN_H
#define ATEAMS_STATISTICS_CHAIN_H

#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <generator>


namespace ATEAMS {
	namespace statistics {
		/**
		 * @class Chain
		 * @brief Simulates the Markov chains specified by each @ref models::Model.
		 * 
		 * @var Chain::model
		 * 	(Pointer to) a Model.
		 * 
		 * @var Chain::state
		 * 	A @ref models::ModelState.
		 * 
		 * @var Chain::options
		 * 	Multithreaded computing options @ref arithmetic::ThreadOptions.
		 * 
		 * @var Chain::steps
		 * 	Number of iterations.
		 */
		template <typename RingLike, template <typename> typename VectorLike>
		class Chain {
			public:
				models::Model<RingLike,VectorLike>* model;
				models::ModelState<RingLike,VectorLike> state;
				arithmetic::ThreadOptions options;

				int steps;

				/**
				 * @brief Constructor; uses default @ref arithmetic::ThreadOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 */
				Chain(models::Model<RingLike,VectorLike>* model, int steps) {
					this->model = model;
					this->steps = steps;

					arithmetic::ThreadOptions options;
					this->options = options;

					models::ModelState<RingLike,VectorLike> state;
					this->state = state;
				};

				/**
				 * @brief Constructor; use user-provided @ref ATEAMS::arithmetic::ThreadOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 * @param options User-provided compute options.
				 */
				Chain(models::Model<RingLike,VectorLike>* model, int steps, arithmetic::ThreadOptions options) {
					this->model = model;
					this->steps = steps;
					this->options = options;

					models::ModelState<RingLike,VectorLike> state;
					this->state = state;
				};

				/**
				 * @brief Model-generic iterator.
				 * 
				 * @code
				 * 	for (models::ModelState<Zp,SparseVector> State : M.simulate()) {
				 * 		<do whatever>
				 * 	}
				 * @endcode
				 * 
				 * @returns A `std::generator`.
				 */
				std::generator<models::ModelState<RingLike,VectorLike>> simulate() {
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

