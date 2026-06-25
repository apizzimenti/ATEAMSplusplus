
#ifndef ATEAMS_MODELS_INVADEDCLUSTER_H
#define ATEAMS_MODELS_INVADEDCLUSTER_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <random>
#include <string>

namespace ATEAMS {
	namespace models {
		/**
		 * @struct InvadedClusterParameters
		 * @brief Parameters relevant to simulating the Swendsen--Wang algorithm.
		 * 
		 * @var InvadedClusterParameters::field
		 * Order \f$p\f$ for the finite field \f$\Z / p \Z\f$.
		 * 
		 * @var InvadedClusterParameters::stoppingFunction
		 * Function that consumes a time-step (integer) parameter \f$t\f$ and returns
		 * an integer in \f$(a,b)\f$, where \f$a, \ b \f$ are between \f$0\f$ and
		 * \f$\rank(H_d(X))\f$ telling us when to sample the next set of spins.
		 */
		struct InvadedClusterParameters : ModelParameters {
			int field = 2;
			std::function<int(int)> stoppingFunction;
		};

		/**
		 * @struct InvadedClusterState
		 * @brief Stores the current state of the Swendsen--Wang algorithm.
		 * 
		 * @var InvadedClusterState::cochain
		 * A (sparse) vector representing the cochain \f$f_t \in C^{d-1}(X; \Z/p\Z)\f$
		 * (i.e. a linear map from \f$(d-1)\f$-cells to \f$\Z/p\Z\f$)
		 * sampled at time \f$t\f$.
		 * 
		 * @var InvadedClusterState::includes
		 * A (sparse) vector with nonzero entries indicating the indices of \f$d\f$-cells
		 * included in the percolation subcomplex sampled at time \f$t\f$.
		 */
		struct InvadedClusterState : ModelState {
			ZpVector cochain;
			std::vector<int> includes;
		};

		/**
		 * @class InvadedCluster
		 * @brief Implements the Swendsen--Wang algorithm for Potts lattice gauge theory
		 * 	(PLGT) and the plaquette random-cluster model (PRCM).
		 * 
		 * @var InvadedCluster::parameters
		 * Model parameters.
		 * 
		 * @var InvadedCluster::state
		 * Model state.
		 * 
		 * @var InvadedCluster::field
		 * Finite field over which computations will be performened.
		 * 
		 * @var InvadedCluster::RNG
		 * Mersenne twister random number generator.
		 * 
		 * @var InvadedCluster::intuniform
		 * Uniform distribution over the integers \f$[0,p)\f$, where \f$p\f$ is
		 * the order of the model's finite field.
		 * 
		 * @var InvadedCluster::kind
		 * Model name.
		 */
		class InvadedCluster: public Model<ZpVector> {
			public:
				/**
				 * @brief Constructor.
				 * 
				 * @param complex (Pointer to) a complex.
				 * @param parameters Model parameters.
				 */
				InvadedCluster(ATEAMS::complexes::Complex* complex, InvadedClusterParameters parameters);

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
				ZpVector sample(int t, ATEAMS::arithmetic::ThreadOptions& options) override;

				/**
				 * @brief Initializes \f$f_0\f$ to uniform random elements of \f$\Z/p\Z\f$.
				 */
				void initialize() override;

				/**
				 * @brief Initializes \f$f_0 = c\f$.
				 */
				void initialize(ZpVector c) override;

				InvadedClusterParameters parameters;
				InvadedClusterState state;
				const Zp field;

				std::string kind = "InvadedCluster";

			private:
				std::mt19937 RNG;
				std::uniform_int_distribution<int> intuniform;
				std::vector<int> filtration;
		};
	}
}

#endif

