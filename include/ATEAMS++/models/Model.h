
#ifndef ATEAMS_MODELS_MODEL_H
#define ATEAMS_MODELS_MODEL_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/complexes/Complex.h"
#include "ATEAMS++/arithmetic/options.h"

#include <string>

namespace ATEAMS::models {
	/**
	 * @struct ModelParameters
	 * @brief Typical parameters required for simulating a Model.
	 * 
	 * @var ModelParameters::coefficients
	 * 	@brief Coefficient ring.
	 * 
	 * @var ModelParameters::dimension
	 * 	Percolation subcomplex dimension. Used by every @ref Model.
	 * 
	 * @var ModelParameters::p
	 * 	Bernoulli trial density. Used by @ref Bernoulli.
	 * 
	 * @var ModelParameters::temperatureFunction
	 * 	Function that specifies the (inverse) temperature parameter at the current
	 * 	time-step. Used by @ref Glauber, @ref SwendsenWang.
	 * 
	 * @var ModelParameters::stoppingFunction
	 * 	Function that specifies the number of giant cycles encountered before
	 * 	re-sampling spins. Used by @ref InvadedCluster, @ref Invasion.
	 * 
	 * @var ModelParameters::DEBUG
	 * Are we debugging the @ref Model?
	 */
	struct ModelParameters {
		// Used by: SwendsenWang, InvadedCluster, Invasion, Glauber
		Ring* coefficients;

		// Used by: Bernoulli
		float p;

		// Used by: SwendsenWang, Glauber
		std::function<double(int)> temperatureFunction;

		// Used by: InvadedCluster, Invasion
		std::function<int(int)> stoppingFunction;

		// Universal.
		int dimension;
		bool DEBUG = false;
	};

	/**
	 * @struct ModelState
	 * @brief Maintains state for this Model.
	 * 
	 * @tparam RingLike Coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @tparam VectorLike Storage unit compatible with @p RingLike, like
	 *  @ref DenseVector or @ref SparseVector.
	 * 
	 * @var ModelState::cochain
	 * 	@brief Current cochain \f$ f_t \f$. Used by @ref InvadedCluster, @ref SwendsenWang.
	 * 
	 * @var ModelState::includes
	 * 	@brief Percolation subcomplex \f$ P_t \f$. Used by @ref Bernoulli, @ref InvadedCluster,
	 * 	@ref Invasion, @ref SwendsenWang.
	 * 
	 * @var ModelState::essential
	 * 	@brief Times at which essential cycles of \f$ P_t \f$ were born. Used by
	 * 	@ref Bernoulli, @ref Invasion.
	 * 
	 * @var ModelState::rank
	 * 	@brief Rank of the \f$d\f$th persistent homology group \f$PH_d(P_t)\f$. Used by
	 * 	@ref Bernoulli.
	 * 
	 * @var ModelState::energy
	 * 	@brief Energy of the current cochain, i.e. \f$\H(f_t)\f$. Used by @ref Glauber.
	 * 
	 * @var ModelState::t
	 * 	@brief Current time-step. Used by every @ref Model.
	 */
	template <typename RingLike, template <typename> typename VectorLike>
	struct ModelState {
		// Used by: SwendsenWang, InvadedCluster,
		VectorLike<RingLike> cochain;

		// Used by: SwendsenWang, InvadedCluster, Invasion, Bernoulli
		std::vector<int> includes;

		// Used by: Bernoulli, Invasion
		std::vector<int> essential;

		// Used by: Bernoulli
		int rank;

		// Used by: Glauber
		int energy;

		// Used by: everyone
		int t;
	};

	/**
	 * @class Model
	 * @brief Abstract (template) class for various Models.
	 * 
	 * @tparam RingLike Coefficient @ref Ring, like @ref Zp or @ref Q.
	 * @tparam VectorLike Storage unit compatible with @p RingLike, like
	 *  @ref DenseVector or @ref SparseVector.
	 * 
	 * @var Model::complex
	 *  @brief Pointer to a @ref ATEAMS::complexes::Complex.
	 * 
	 * @var Model::coefficients
	 *  @brief Coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @var Model::dimension
	 * 	@brief Subcomplex dimension.
	 * 
	 * @var Model::DEBUG
	 * 	@brief Are we debugging this Model? (Default is `false`.)
	 */
	template <typename RingLike, template <typename> typename VectorLike>
	class Model {
		public:
			ATEAMS::complexes::Complex<RingLike>* complex;

			Ring* coefficients;
			int dimension;
			bool DEBUG = false;

			/** Constructor. */
			Model(
				Ring* R,
				int dimension,
				bool DEBUG=false

			) : coefficients(R), dimension(dimension), DEBUG(DEBUG) { };

			/**
			 * @brief Draw a sample from the model.
			 * 
			 * @param t When simulated by a Markov chain, \f$t\f$ is the current
			 * 	time-step.
			 * @param state Keeps track of state.
			 * @param options Options for the multithreaded computing environment.
			 * 
			 * @returns Modified state.
			 */
			virtual ModelState<RingLike,VectorLike> sample(int t, ModelState<RingLike,VectorLike>& state, ATEAMS::arithmetic::ThreadOptions& options) = 0;

			/**
			 * @brief Initializes the state as determined by the model.
			 * 
			 * @param state Keeps track of state.
			 * 
			 * @returns Modified state.
			 */
			virtual ModelState<RingLike,VectorLike> initialize(ModelState<RingLike,VectorLike>& state) = 0;

			/**
			 * @brief Initializes the state as determined by the user.
			 * 
			 * @param c Initial cochain.
			 * @param state Keeps track of state.
			 * 
			 * @returns Modified state.
			 */
			virtual ModelState<RingLike,VectorLike> initialize(VectorLike<RingLike> c, ModelState<RingLike,VectorLike>& state) = 0;
	};
}

#endif
