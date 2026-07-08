
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
		 * @class InvadedCluster
		 * @brief Implements the Swendsen--Wang algorithm for Potts lattice gauge theory
		 * 	(PLGT) and the plaquette random-cluster model (PRCM).
		 * 
		 * @var InvadedCluster::parameters
		 * Model parameters.
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
		template <typename T>
		class InvadedCluster: public Model<T,SparseVector> {
			public:
				/**
				 * @brief Constructor.
				 * 
				 * @param complex (Pointer to) a complex.
				 * @param parameters Model parameters.
				 */
				InvadedCluster(ATEAMS::complexes::Complex<T>* complex, ModelParameters parameters);

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
				ModelState<T,SparseVector> sample(int t, ModelState<T,SparseVector>& state, ATEAMS::arithmetic::ThreadOptions& options) override;

				/**
				 * @brief Initializes \f$f_0\f$ to uniform random elements of \f$\Z/p\Z\f$.
				 */
				ModelState<T,SparseVector> initialize(ModelState<T,SparseVector>& state) override;

				/**
				 * @brief Initializes \f$f_0 = c\f$.
				 */
				ModelState<T,SparseVector> initialize(SparseVector<T> c, ModelState<T,SparseVector>& state) override;

				ModelParameters parameters;
				const Field field;
				std::string kind = "InvadedCluster";

			private:
				std::mt19937 RNG;
				std::uniform_int_distribution<int> intuniform;
				std::vector<int> filtration;
		};
	}
}

#include "ATEAMS++/models/InvadedCluster.tpp"

#endif

