
#ifndef ATEAMS_STATISTICS_CHAIN_H
#define ATEAMS_STATISTICS_CHAIN_H

#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <generator>


namespace ATEAMS {
	namespace statistics {

		template <typename ModelType>
		class Chain {
			public:
				ModelType* model;
				ATEAMS::models::ModelState* state;
				ATEAMS::arithmetic::ThreadOptions options;

				int steps;

				Chain(ModelType* model, int steps) {
					this->model = model;
					this->steps = steps;

					ATEAMS::arithmetic::ThreadOptions options;
					this->options = options;
				};

				Chain(ModelType* model, int steps, ATEAMS::arithmetic::ThreadOptions options) {
					this->model = model;
					this->steps = steps;
					this->options = options;
				};

				template <typename State>
				std::generator<State*> run() {
					// Spin up multithreading environment.
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

