
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
		 * 	A @ref ATEAMS::models::ModelState.
		 * 
		 * @var Chain::options
		 * 	Multithreaded computing options @ref ATEAMS::arithmetic::ThreadOptions.
		 * 
		 * @var Chain::steps
		 * 	Number of iterations.
		 */
		template <typename T=ATEAMS::ff, template <typename> typename ContainerType=ATEAMS::SparseVector>
		class Chain {
			public:
				ATEAMS::models::Model<T,ContainerType>* model;
				ATEAMS::models::ModelState<T,ContainerType> state;
				ATEAMS::arithmetic::ThreadOptions options;

				int steps;

				/**
				 * @brief Constructor; uses default @ref ATEAMS::arithmetic::ThreadOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 */
				Chain(ATEAMS::models::Model<T,ContainerType>* model, int steps) {
					this->model = model;
					this->steps = steps;

					ATEAMS::arithmetic::ThreadOptions options;
					this->options = options;

					ATEAMS::models::ModelState<T,ContainerType> state;
					this->state = state;
				};

				/**
				 * @brief Constructor; use user-provided @ref ATEAMS::arithmetic::ThreadOptions.
				 * 
				 * @param model (Pointer to) a Model instance.
				 * @param steps Number of iterations.
				 * @param options User-provided compute options.
				 */
				Chain(ATEAMS::models::Model<T,ContainerType>* model, int steps, ATEAMS::arithmetic::ThreadOptions options) {
					this->model = model;
					this->steps = steps;
					this->options = options;

					ATEAMS::models::ModelState<T,ContainerType> state;
					this->state = state;
				};

				/**
				 * @brief Model-generic iterator.
				 * 
				 * @code
				 * 	for (models::ModelState<ATEAMS::ff,ATEAMS::SparseVector> State : M.simulate()) {
				 * 		<do whatever>
				 * 	}
				 * @endcode
				 * 
				 * @returns A `std::generator`.
				 */
				std::generator<ATEAMS::models::ModelState<T,ContainerType>> simulate() {
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

