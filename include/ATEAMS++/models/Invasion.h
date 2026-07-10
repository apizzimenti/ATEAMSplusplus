
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
		 * @class Invasion
		 * @brief Implements the invaded-cluster algorithm for Potts lattice gauge theory
		 * 	(PLGT) and the plaquette random-cluster model (PRCM).
		 * @tparam RingLike A coefficient @ref Ring type, like @ref Zp or @ref Q.
		 * 
		 * **Example Usage**
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
		 * models::Invasion<CC::RingType> invasion(&plex, params);
		 * @endcode
		 * 
		 * @var Invasion::stoppingFunction
		 * 	@brief Specifies the number of homological percolation events before
		 * 	re-sampling the cochain. See @ref ATEAMS::statistics::stopInvadingAt.
		 * 
		 * @var Invasion::name
		 * 	@brief Human-readable name.
		 */
		template <typename RingLike>
		class Invasion: public Model<RingLike,DenseVector> {
			public:
				std::function<int(int)> stoppingFunction;
				std::string name = "Invasion";

				/**
				 * @brief Exposed coefficient ring type. See @ref Ring.
				 */
				using RingType = RingLike;

				/**
				 * @brief Exposed vector storage type. See @ref DenseVector.
				 */
				template <typename R>
				using VectorType = DenseVector<R>;

				/**
				 * @brief Constructor.
				 * 
				 * @param complex (Pointer to) a complex.
				 * @param parameters Model parameters.
				 */
				Invasion(ATEAMS::complexes::Complex<RingLike>* complex, ModelParameters parameters);

				/**
				 * @brief Constructor.
				 * 
				 * @param complex (Pointer to) a complex.
				 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
				 * @param dimension Subcomplex dimension.
				 * @param stoppingFunction A function that consumes a time-step and
				 * 	returns an integer number of giant cycles to be observed before
				 * 	invasion percolation is halted.
				 * @param DEBUG _(Optional, default `false`)_ Are we debugging this Model?
				 */
				Invasion(
					ATEAMS::complexes::Complex<RingLike>* complex,
					Ring* R,
					int dimension,
					std::function<int(int)> stoppingFunction,
					bool DEBUG=false
				);

				/**
				 * @brief Implements the plaquette invasion-percolation algorithm,
				 * which creates a random \f$d\f$-subcomplex \f$X\f$ by (uniformly 
				 * randomly) inserting \f$d\f$-cells into \f$X\f$ until \f$\rank(H_d(X))\f$
				 * is "large enough" (set by the user).
				 * 
				 * @param t Time step.
				 * @param state Model state.
				 * @param options Multithreaded computing environment options.
				 * 
				 * @return Model state with @ref ModelState::includes, @ref ModelState::essential,
				 * and @ref ModelState::t updated.
				 */
				ModelState<RingLike,DenseVector> sample(
					int t,
					ModelState<RingLike,DenseVector>& state,
					ATEAMS::arithmetic::ThreadOptions& options
				) override;

				/** @brief Initialization; superfluous. */
				ModelState<RingLike,DenseVector> initialize(ModelState<RingLike,DenseVector>& state) override { return state; };

				/** @brief Initialization; superfluous. */
				ModelState<RingLike,DenseVector> initialize(std::vector<RingLike> c, ModelState<RingLike,DenseVector>& state) override { return state; };

			private:
				std::mt19937 RNG;
				std::uniform_int_distribution<int> intuniform;
				std::vector<int> indices;
				std::vector<int> filtration;
		};
	}
}

#include "ATEAMS++/models/Invasion.tpp"

#endif

