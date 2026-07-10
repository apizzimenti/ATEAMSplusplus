
#ifndef ATEAMS_MODELS_INVADEDCLUSTER_H
#define ATEAMS_MODELS_INVADEDCLUSTER_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <random>
#include <string>

namespace ATEAMS::models {
	/**
	 * @class InvadedCluster
	 * @brief Implements the invaded-cluster algorithm.
	 * @tparam RingLike A coefficient @ref Ring type, like @ref Zp or @ref Q.
	 * 
	 * @code
	 * // Complexes must be constructed with the same underlying coefficient
	 * // ring as the model. If the ring is known (e.g. for Bernoulli percolation it
	 * // is always ATEAMS::Z2), you can name it explicitly; otherwise, the
	 * // model exposes the ring type after it is declared. You can also create
	 * // the complex *first* using a specific ring, which is exposed under the
	 * // same type alias (i.e. ::RingType).
	 * using CC = complexes::Cubical<Zp>;
	 * 
	 * CC plex({4,4,4,4});
	 * CC::RingType RR(3); // equivalent to Zp RR(3);
	 * 
	 * models::ModelParameters params;
	 * params.dimension = 2;
	 * params.coefficients = &RR;
	 * params.stoppingFunction = statistics::stopInvadingAt({3,4});
	 * 
	 * models::InvadedCluster<CC::RingType> IC(&plex, params);
	 * @endcode
	 * 
	 * @var InvadedCluster::stoppingFunction
	 * 	@brief Specifies the number of homological percolation events before
	 * 	re-sampling the cochain. See @ref ATEAMS::statistics::stopInvadingAt.
	 * 
	 * @var InvadedCluster::name
	 * 	@brief Human-readable name.
	 */
	template <typename RingLike>
	class InvadedCluster: public Model<RingLike,SparseVector> {
		public:
			std::function<int(int)> stoppingFunction;
			std::string name = "Invaded-Cluster";

			/**
			 * @brief Exposed coefficient ring type. See @ref Ring.
			 */
			using RingType = RingLike;

			/**
			 * @brief Exposed vector storage type. See @ref SparseVector.
			 */
			template <typename R>
			using VectorType = SparseVector<R>;

			/**
			 * @brief Constructor.
			 * 
			 * @param complex (Pointer to) a complex.
			 * @param parameters Model parameters.
			 */
			InvadedCluster(
				ATEAMS::complexes::Complex<RingLike>* complex,
				ModelParameters parameters
			);

			/**
			 * @brief Constructor.
			 * 
			 * Using this constructor, the `stoppingFunction` is automatically
			 * 	to stop invading after half (or half + 1) of the giant cycles have
			 * 	been found; see @ref ATEAMS::statistics::stopInvadingAt.
			 * 
			 * @param complex (Pointer to) a complex.
			 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
			 * @param dimension Subcomplex dimension.
			 * @param stoppingFunction A function that consumes a time-step and
			 * 	returns an integer number of giant cycles to be observed before
			 * 	invasion percolation is halted.
			 * @param DEBUG _(Optional, default `false`)_ Are we debugging this Model?
			 */
			InvadedCluster(
				ATEAMS::complexes::Complex<RingLike>* complex,
				Ring* R,
				int dimension,
				std::function<int(int)> stoppingFunction,
				bool DEBUG=false
			);

			/**
			 * @brief Implements the plaquette invaded-cluster algorithm,
			 * 	which inserts satisfied cells (i.e. cells \f$x\f$ such that
			 * 	\f$(\delta^{d-1}f_t)(x) = 0\f$) in a uniform random order until
			 * 	some number of giant cycles are created, forming a subcomplex
			 * 	\f$P_t\f$. Then samples \f$ f_{t+1} \sim \PK(- \mid P_t)\f$.
			 * 
			 * @param t Time step.
			 * @param state Model state.
			 * @param options Multithreaded computing environment options.
			 * 
			 * @return Model state with @ref ModelState::cochain, @ref ModelState::includes,
			 * @ref ModelState::essential, and @ref ModelState::t updated.
			 */
			ModelState<RingLike,SparseVector> sample(
				int t,
				ModelState<RingLike,SparseVector>& state,
				ATEAMS::arithmetic::ThreadOptions& options
			) override;

			/**
			 * @brief Initializes \f$f_0\f$ to uniform random element of \f$\Z/p\Z\f$
			 * or \f$ \Q \f$.
			 * 
			 * @param state Model state.
			 * 
			 * @returns Model state.
			 */
			ModelState<RingLike,SparseVector> initialize(ModelState<RingLike,SparseVector>& state) override;

			/**
			 * @brief Initializes \f$f_0 = c\f$.
			 * 
			 * @param c Cochain to which \f$f_0\f$ is initialized.
			 * @param state Model state.
			 * 
			 * @returns Model state.
			 */
			ModelState<RingLike,SparseVector> initialize(
				SparseVector<RingLike> c,
				ModelState<RingLike,SparseVector>& state
			) override;

		private:
			std::mt19937 RNG;
			std::uniform_int_distribution<int> intuniform;
			std::vector<int> filtration;
	};
}

#include "ATEAMS++/models/InvadedCluster.tpp"

#endif

