
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
		 */
		class Invasion: public Model<ff,DenseVector> {
			public:
				/** Coefficients are of type @ref ATEAMS::ff. */
				typedef ff dt;

				/** Vectors are of type @ref ATEAMS::DenseVector. */
				template <typename R>
				using st = DenseVector<R>;

				/**
				 * @brief Constructor.
				 * 
				 * @param complex (Pointer to) a complex.
				 * @param parameters Model parameters.
				 */
				Invasion(ATEAMS::complexes::Complex<ff>* complex, ModelParameters parameters);

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
				ModelState<ff,DenseVector> sample(int t, ModelState<ff,DenseVector>& state, ATEAMS::arithmetic::ThreadOptions& options) override;

				/** @brief Initialization; superfluous. */
				ModelState<ff,DenseVector> initialize(ModelState<ff,DenseVector>& state) override { return state; };

				/** @brief Initialization; superfluous. */
				ModelState<ff,DenseVector> initialize(std::vector<ff> c, ModelState<ff,DenseVector>& state) override { return state; };

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

