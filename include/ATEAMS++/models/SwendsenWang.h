
#ifndef ATEAMS_MODELS_SWENDSENWANG_H
#define ATEAMS_MODELS_SWENDSENWANG_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <random>
#include <string>

namespace ATEAMS::models {
	/**
	 * @class SwendsenWang
	 * @brief Implements the Swendsen--Wang algorithm for Potts lattice gauge theory
	 * 	(PLGT) and the plaquette random-cluster model (PRCM).
	 * 
	 * @tparam RingLike A coefficient @ref Ring type, like @ref Zp or @ref Q.
	 * 
	 * **Example Usage**
	 * 
	 * @code
	 * // Complexes must be constructed with the same underlying coefficient
	 * // ring as the model. If the ring is known (e.g. for Bernoulli percolation, the ring is always ATEAMS::Z2) you can name it explicitly; otherwise, the
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
	 * params.stoppingFunction = statistics::selfdual(&RR);
	 * 
	 * models::SwendsenWang<CC::RingType> SW(&plex, params);
	 * @endcode
	 * 
	 * @var SwendsenWang::temperatureFunction
	 * 	@brief Specifies an inverse temperature at each time-step; see @ref ATEAMS::statistics::selfdual.
	 * 
	 * @var SwendsenWang::name
	 * 	@brief Human-readable name.
	 */
	template <typename RingLike>
	class SwendsenWang: public Model<RingLike,SparseVector> {
		public:
			std::function<double(int)> temperatureFunction;
			std::string name = "Swendsen--Wang";

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
			 * Using this constructor, typical model parameters can be passed
			 * as a @ref ATEAMS::models::ModelParameters object.
			 * 
			 * @param complex (Pointer to) a complex.
			 * @param parameters Model parameters.
			 */
			SwendsenWang(
				ATEAMS::complexes::Complex<RingLike>* complex,
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
			SwendsenWang(
				ATEAMS::complexes::Complex<RingLike>* complex,
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
			SwendsenWang(
				ATEAMS::complexes::Complex<RingLike>* complex,
				Ring* R,
				int dimension,
				std::function<double(int)> temperatureFunction,
				bool DEBUG=false
			);

			/**
			 * @brief Model sampler.
			 * 
			 * Alternates between sampling the conditionals \f$P_t \sim \PK(- \mid f_t)\f$
			 * 	and \f$f_{t+1} \sim \PK(- \mid P_{t+1})\f$ of the coupling \f$\PK\f$ between
			 * 	the PLGT and PRCM.
			 * 
			 * @param t Time step.
			 * @param state State.
			 * @param options Multithreaded computing environment options.
			 * 
			 * @return Model state, with @ref ModelState::cochain, @ref ModelState::includes,
			 * and @ref ModelState::t updated.
			 */
			ModelState<RingLike,SparseVector> sample(int t, ModelState<RingLike,SparseVector>& state, ATEAMS::arithmetic::ThreadOptions& options) override;

			/**
			 * @brief Initializes \f$f_0\f$ to uniform random element of \f$\Z/p\Z\f$
			 * or \f$ \Q \f$.
			 * 
			 * @param state State.
			 * 
			 * @returns Modified State.
			 */
			ModelState<RingLike,SparseVector> initialize(ModelState<RingLike,SparseVector>& state) override;

			/**
			 * @brief Initializes \f$f_0 = c\f$.
			 * 
			 * @param c Cochain to which \f$f_0\f$ is initialized.
			 * @param state State.
			 * 
			 * @returns Modified State.
			 */
			ModelState<RingLike,SparseVector> initialize(SparseVector<RingLike> c, ModelState<RingLike,SparseVector>& state) override;

		private:
			std::mt19937 RNG;
			std::uniform_real_distribution<double> unituniform;
			std::uniform_int_distribution<int> intuniform;
			std::set<int> include;
	};
}

#include "ATEAMS++/models/SwendsenWang.tpp"

#endif

