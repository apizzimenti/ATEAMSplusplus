
#ifndef ATEAMS_MODELS_BERNOULLI_H
#define ATEAMS_MODELS_BERNOULLI_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/persistence.h"

#include <random>
#include <set>

namespace ATEAMS {
	namespace models {
		/**
		 * @class Bernoulli
		 * @brief Implements Bernoulli percolation.
		 */
		class Bernoulli : public Model<FINITE,DenseVector> {
			public:
				/** Coefficients are of type @ref ATEAMS::FINITE. */
				typedef FINITE dt;

				/** Vectors are of type @ref ATEAMS::DenseVector. */
				template <typename R>
				using st = DenseVector<R>;

				/**
				 * @brief Constructor.
				 * @param complex (Pointer to) a Complex.
				 * @param parameters Model parameters.
				 */
				Bernoulli(ATEAMS::complexes::Complex<FINITE>* complex, ModelParameters parameters);

				/**
				 * @brief Directly samples Bernoulli percolation at the desired dimension.
				 * @param t Time step.
				 * @param state Tracks the state of the model. See @ref ATEAMS::statistics::Chain.
				 * @param options Multithreaded computing environment options.
				 * 
				 * @returns @ref ModelState with @ref ModelState::essential,
				 * @ref ModelState::includes, and @ref ModelState::t modified.
				 */
				ModelState<FINITE,DenseVector> sample(int t, ModelState<FINITE,DenseVector>& state, ATEAMS::arithmetic::ThreadOptions& options) override;

				/** @brief Initialization; superfluous. */
				ModelState<FINITE,DenseVector> initialize(ModelState<FINITE,DenseVector>& state) override { return state; };

				/** @brief Initialization; superfluous. */
				ModelState<FINITE,DenseVector> initialize(std::vector<FINITE> c, ModelState<FINITE,DenseVector>& state) override { return state; };
			private:
				std::mt19937 RNG;
				std::uniform_real_distribution<double> unituniform;
				const Field field;
				
				std::vector<int> filtration;
		};
	}
}

#include "ATEAMS++/models/Bernoulli.tpp"

#endif
