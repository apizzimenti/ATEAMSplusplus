
#ifndef ATEAMS_MODELS_SWENDSENWANG_H
#define ATEAMS_MODELS_SWENDSENWANG_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <random>

namespace ATEAMS {
	namespace models {
		/**
		 * @struct SwendsenWangParameters
		 * @brief Parameters relevant to simulating the Swendsen--Wang algorithm.
		 * 
		 * @var SwendsenWangParameters::field
		 * Order \f$p\f$ for the finite field \f$\Z / p \Z\f$.
		 * 
		 * @var SwendsenWangParameters::temperatureFunction
		 * Function that consumes a time-step (integer) parameter \f$t\f$, and returns
		 * a temperature parameter \f$\beta \in [0, \infty)\f$ that is used to compute
		 * the probability \f$p = 1-e^{-\beta}\f$ of including any individual
		 * \f$(d-1)\f$-cell in the sampled subcomplex.
		 */
		struct SwendsenWangParameters : ModelParameters {
			int field = 2;
			std::function<double(int)> temperatureFunction;
		};

		/**
		 * @struct SwendsenWangState
		 * @brief Stores the current state of the Swendsen--Wang algorithm.
		 * 
		 * @var SwendsenWangState::cochain
		 * A (sparse) vector representing the cochain \f$f_t \in C^{d-1}(X; \Z/p\Z)\f$
		 * (i.e. a linear map from \f$(d-1)\f$-cells to \f$\Z/p\Z\f$)
		 * sampled at time \f$t\f$.
		 * 
		 * @var SwendsenWangState::includes
		 * A (sparse) vector with nonzero entries indicating the indices of \f$d\f$-cells
		 * included in the percolation subcomplex sampled at time \f$t\f$.
		 */
		struct SwendsenWangState : ModelState {
			ZpVector cochain;
			ZpVector includes;
		};

		/**
		 * @class SwendsenWang
		 * @brief Implements the Swendsen--Wang algorithm for Potts lattice gauge theory
		 * 	(PLGT) and the plaquette random-cluster model (PRCM).
		 * 
		 * @var SwendsenWang::parameters
		 * Model parameters.
		 * 
		 * @var SwendsenWang::state
		 * Model state.
		 * 
		 * @var SwendsenWang::field
		 * Finite field over which computations will be performened.
		 * 
		 * @var SwendsenWang::RNG
		 * Mersenne twister random number generator.
		 * 
		 * @var SwendsenWang::unituniform
		 * Uniform distribution over \f$(0,1)\f$.
		 * 
		 * @var SwendsenWang::intuniform
		 * Uniform distribution over the integers \f$[0,p)\f$, where \f$p\f$ is
		 * the order of the model's finite field.
		 */
		class SwendsenWang: public Model {
			public:
				SwendsenWang(ATEAMS::complexes::Complex* complex, SwendsenWangParameters parameters);

				ZpVector sample(int t, ATEAMS::arithmetic::ThreadOptions& options) override;

				void initialize() override;
				void initialize(ZpVector c) override;

				SwendsenWangParameters parameters;
				SwendsenWangState state;
				const Zp field;

			private:
				std::mt19937 RNG;
				std::uniform_real_distribution<double> unituniform;
				std::uniform_int_distribution<int> intuniform;
		};
	}
}

#endif

