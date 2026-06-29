
#ifndef ATEAMS_MODELS_INVASION_H
#define ATEAMS_MODELS_INVASION_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <random>
#include <string>

namespace ATEAMS {
	namespace models {
		/**
		 * @struct InvasionParameters
		 * @brief Parameters relevant to simulating the invaded-cluster algorithm.
		 * 
		 * @var InvasionParameters::field
		 * Order \f$p\f$ for the finite field \f$\Z / p \Z\f$.
		 * 
		 * @var InvasionParameters::stoppingFunction
		 * Function that consumes a time-step (integer) parameter \f$t\f$ and returns
		 * an integer in \f$(a,b)\f$, where \f$a, \ b \f$ are between \f$0\f$ and
		 * \f$\rank(H_d(X))\f$ telling us when to sample the next set of spins.
		 */
		struct InvasionParameters : ModelParameters {
			int field = 2;
			std::function<int(int)> stoppingFunction;
		};

		/**
		 * @struct InvasionState
		 * @brief Stores the current state of the invaded-cluster algorithm.
		 * 
		 * @var InvasionState::includes
		 * A vector with nonzero entries indicating the indices of \f$d\f$-cells
		 * included in the percolation subcomplex sampled at time \f$t\f$, in
		 * insertion order.
		 * 
		 * @var InvasionState::essential
		 * A vector with the times at which giant cycles were encountered.
		 */
		struct InvasionState : ModelState {
			std::vector<int> includes;
			std::vector<int> essential;
		};

		/**
		 * @class Invasion
		 * @brief Implements the invaded-cluster algorithm for Potts lattice gauge theory
		 * 	(PLGT) and the plaquette random-cluster model (PRCM).
		 * 
		 * @var Invasion::parameters
		 * Model parameters.
		 * 
		 * @var Invasion::state
		 * Model state.
		 * 
		 * @var Invasion::field
		 * Finite field over which computations will be performened.
		 * 
		 * @var Invasion::RNG
		 * Mersenne twister random number generator.
		 * 
		 * @var Invasion::intuniform
		 * Uniform distribution over the integers \f$[0,p)\f$, where \f$p\f$ is
		 * the order of the model's finite field.
		 * 
		 * @var Invasion::kind
		 * Model name.
		 */
		class Invasion: public Model<std::vector<int>> {
			public:
				/**
				 * @brief Constructor.
				 * 
				 * @param complex (Pointer to) a complex.
				 * @param parameters Model parameters.
				 */
				Invasion(ATEAMS::complexes::Complex* complex, InvasionParameters parameters);

				/**
				 * @brief Implements the plaquette invasion-percolation algorithm,
				 * which creates a random \f$d\f$-subcomplex \f$X\f$ by (uniformly 
				 * randomly) inserting \f$d\f$-cells into \f$X\f$ until \f$\rank(H_d(X))\f$
				 * is "large enough" (set by the user).
				 * 
				 * @param t Time step.
				 * @param options Multithreaded computing environment options.
				 * 
				 * @return A list of \f$d\f$-cell indices in insertion order, up
				 * to the time at which \f$\rank(H_d(X))\f$ exceeds the user-specified
				 * cutoff.
				 */
				std::vector<int> sample(int t, ATEAMS::arithmetic::ThreadOptions& options) override;

				/** @brief Initialization; superfluous. */
				void initialize() override { };

				/** @brief Initialization; superfluous. */
				void initialize(std::vector<int> c) override { };

				InvasionParameters parameters;
				InvasionState state;
				const Zp field;

				std::string kind = "Invasion";

			private:
				std::mt19937 RNG;
				std::uniform_int_distribution<int> intuniform;
				std::vector<int> indices;
				std::vector<int> filtration;
		};
	}
}

#endif

