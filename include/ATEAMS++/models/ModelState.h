
#ifndef ATEAMS_MODELS_MODELSTATE_H
#define ATEAMS_MODELS_MODELSTATE_H

#include "ATEAMS++/common.h"

namespace ATEAMS::models {
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
}

#endif
