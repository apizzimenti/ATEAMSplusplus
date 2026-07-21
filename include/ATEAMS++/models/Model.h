
#ifndef ATEAMS_MODELS_MODEL_H
#define ATEAMS_MODELS_MODEL_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/ModelState.h"
#include "ATEAMS++/models/ModelParameters.h"
#include "ATEAMS++/complexes/Complex.h"
#include "ATEAMS++/arithmetic/options.h"

#include <string>

namespace ATEAMS::models {
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
	 *  @brief Pointer to a coefficient @ref Ring, like @ref Zp or @ref Q.
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

			template <typename R>
			using VectorType = VectorLike<R>;
			using RingType = RingLike;
			using State = ModelState<RingLike,VectorLike>;

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
			virtual State sample(int t, State& state, ATEAMS::arithmetic::ComputeOptions<RingLike>& options) = 0;

			/**
			 * @brief Initializes the state as determined by the model.
			 * 
			 * @param state Keeps track of state.
			 * 
			 * @returns Modified state.
			 */
			virtual State initialize(State& state) = 0;

			/**
			 * @brief Initializes the state as determined by the user.
			 * 
			 * @param c Initial cochain.
			 * @param state Keeps track of state.
			 * 
			 * @returns Modified state.
			 */
			virtual State initialize(VectorLike<RingLike> c, State& state) = 0;
	};
}

#endif
