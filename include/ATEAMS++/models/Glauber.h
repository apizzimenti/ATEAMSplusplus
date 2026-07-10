
#ifndef ATEAMS_MODELS_GLAUBER_H
#define ATEAMS_MODELS_GLAUBER_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <random>
#include <string>

namespace ATEAMS::models {
	/**
	 * @class Glauber
	 * @brief Implements the Glauber dynamics algorithm.
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
	 * params.temperatureFunction = statistics::selfdual(&RR);
	 * 
	 * models::Glauber<CC::RingType> glauber(&plex, params);
	 * @endcode
	 * 
	 * @var Glauber::temperatureFunction
	 * 	@brief Specifies an inverse temperature at each time-step; see @ref ATEAMS::statistics::selfdual.
	 * 
	 * @var Glauber::name
	 * 	@brief Human-readable name.
	 */
	template <typename RingLike>
	class Glauber: public Model<RingLike,SparseVector> {
		public:
			std::function<double(int)> temperatureFunction;
			std::string name = "Glauber";

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
			Glauber(
				complexes::Complex<RingLike>* complex,
				ModelParameters parameters
			);

			/**
			 * @brief Constructor.
			 * 
			 * Using this constructor, the `temperatureFunction` is automatically
			 * 	set to the inverse temperature parameter admitting the self-dual
			 * 	point according to the modulus of the @ref Ring; see
			 *  @ref ATEAMS::statistics::selfdual.
			 * 
			 * @param complex (Pointer to) a complex.
			 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
			 * @param dimension Subcomplex dimension.
			 * @param DEBUG _(Optional, default `false`)_ Are we debugging this Model?
			 */
			Glauber(
				complexes::Complex<RingLike>* complex,
				Ring* R,
				int dimension,
				bool DEBUG=false
			);

			/**
			 * @brief Constructor.
			 * 
			 * @param complex (Pointer to) a complex.
			 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
			 * @param dimension Subcomplex dimension.
			 * @param temperatureFunction A function that consumes an integer time \f$t\f$ and
			 * 	returns a double, representing the _inverse temperature_ at time \f$t\f$.
			 * @param DEBUG _(Optional, default `false`)_ Are we debugging this Model?
			 */
			Glauber(
				complexes::Complex<RingLike>* complex,
				Ring* R,
				int dimension,
				std::function<double(int)> temperatureFunction,
				bool DEBUG=false
			);

			/**
			 * @brief Switches the spin of a uniform random cell.
			 * 
			 * @param t Time step.
			 * @param state Tracks state; see @ref ATEAMS::statistics::Chain.
			 * @param options Multithreaded computing environment options.
			 * 
			 * @return State with modified @ref ModelState::cochain,
			 * @ref ModelState::energy, and @ref ModelState::t.
			 */
			ModelState<RingLike,SparseVector> sample(
				int t,
				ModelState<RingLike,SparseVector>& state,
				arithmetic::ThreadOptions& options
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
			ModelState<RingLike,SparseVector> initialize(SparseVector<RingLike> c, ModelState<RingLike,SparseVector>& state) override;

		private:
			std::mt19937 RNG;
			std::uniform_real_distribution<double> unituniform;
			std::uniform_int_distribution<int> intuniform;
			std::uniform_int_distribution<int> indexuniform;
	};
}

#include "ATEAMS++/models/Glauber.tpp"

#endif

