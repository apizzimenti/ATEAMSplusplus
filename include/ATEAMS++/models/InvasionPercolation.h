
#ifndef ATEAMS_MODELS_INVASIONPERCOLATION_H
#define ATEAMS_MODELS_INVASIONPERCOLATION_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <random>
#include <string>

namespace ATEAMS {
	namespace models {
		/**
		 * @struct InvasionPercolationParameters
		 * @brief Parameters relevant to simulating the invaded-cluster algorithm.
		 * 
		 * @var InvasionPercolationParameters::field
		 * Order \f$p\f$ for the finite field \f$\Z / p \Z\f$.
		 * 
		 * @var InvasionPercolationParameters::stoppingFunction
		 * Function that consumes a time-step (integer) parameter \f$t\f$ and returns
		 * an integer in \f$(a,b)\f$, where \f$a, \ b \f$ are between \f$0\f$ and
		 * \f$\rank(H_d(X))\f$ telling us when to sample the next set of spins.
		 */
		struct InvasionPercolationParameters : ModelParameters {
			int field = 2;
			std::function<int(int)> stoppingFunction;
		};

		/**
		 * @struct InvasionPercolationState
		 * @brief Stores the current state of the invaded-cluster algorithm.
		 * 
		 * @var InvasionPercolationState::includes
		 * A (sparse) vector with nonzero entries indicating the indices of \f$d\f$-cells
		 * included in the percolation subcomplex sampled at time \f$t\f$.
		 */
		struct InvasionPercolationState : ModelState {
			std::vector<int> includes;
			std::vector<int> essential;
		};

		/**
		 * @class InvasionPercolation
		 * @brief Implements the invaded-cluster algorithm for Potts lattice gauge theory
		 * 	(PLGT) and the plaquette random-cluster model (PRCM).
		 * 
		 * @var InvasionPercolation::parameters
		 * Model parameters.
		 * 
		 * @var InvasionPercolation::state
		 * Model state.
		 * 
		 * @var InvasionPercolation::field
		 * Finite field over which computations will be performened.
		 * 
		 * @var InvasionPercolation::RNG
		 * Mersenne twister random number generator.
		 * 
		 * @var InvasionPercolation::intuniform
		 * Uniform distribution over the integers \f$[0,p)\f$, where \f$p\f$ is
		 * the order of the model's finite field.
		 * 
		 * @var InvasionPercolation::kind
		 * Model name.
		 */
		class InvasionPercolation: public Model<std::vector<int>> {
			public:
				/**
				 * @brief Constructor.
				 * 
				 * @param complex (Pointer to) a complex.
				 * @param parameters Model parameters.
				 */
				InvasionPercolation(ATEAMS::complexes::Complex* complex, InvasionPercolationParameters parameters);

				/**
				 * @brief Implements the plaquette invaded-cluster algorithm,
				 * 	which inserts satisfied cells (i.e. cells \f$x\f$ such that
				 * 	\f$(\delta^{d-1}f_t)(x) = 0\f$) in a uniform random order until
				 * 	some number of giant cycles are created, forming a subcomplex
				 * 	\f$P_t\f$. Then samples \f$ f_{t+1} \sim \PK(- \mid P_t)\f$.
				 * 
				 * @param t Time step.
				 * @param options Multithreaded computing environment options.
				 * 
				 * @return The sample \f$f_{t+1}\f$.
				 */
				std::vector<int> sample(int t, ATEAMS::arithmetic::ThreadOptions& options) override;

				/** @brief Initialization; superfluous. */
				void initialize() override { };

				/** @brief Initialization; superfluous. */
				void initialize(std::vector<int> c) override { };

				InvasionPercolationParameters parameters;
				InvasionPercolationState state;
				const Zp field;

				std::string kind = "InvasionPercolation";

			private:
				std::mt19937 RNG;
				std::uniform_int_distribution<int> intuniform;
				std::vector<int> indices;
				std::vector<int> filtration;
		};
	}
}

#endif

