
#ifndef ATEAMS_MODELS_BERNOULLI_H
#define ATEAMS_MODELS_BERNOULLI_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/persistence.h"

#include <random>
#include <set>

namespace ATEAMS::models {
	/**
	 * @class Bernoulli
	 * @brief Implements Bernoulli percolation.
	 * 
	 * @var Bernoulli::p
	 * 	@brief Bernoulli trial density \f$p \in [0,1]\f$.
	 * 
	 * @var Bernoulli::name
	 * 	@brief Human-readable name.
	 */
	class Bernoulli : public Model<Z2,DenseVector> {
		public:
			float p;
			std::string name = "Bernoulli";

			/**
			 * @brief Exposed coefficient ring type. See @ref Z2.
			 */
			using CoefficientType = Z2;

			/**
			 * @brief Exposed vector storage type. See @ref SparseVector.
			 */
			template <typename R>
			using VectorType = DenseVector<R>;

			/**
			 * @brief Constructor.
			 * 
			 * Using this constructor, typical model parameters can be passed
			 * as a @ref ATEAMS::models::ModelParameters object.
			 * 
			 * @param complex (Pointer to) a complex.
			 * @param parameters Model parameters.
			 */
			Bernoulli(
				ATEAMS::complexes::Complex<Z2>* complex,
				ModelParameters parameters
			);

			/**
			 * @brief Constructor.
			 * 
			 * Using this constructor, the subcomplex dimension and percolation
			 * density must be passed manually.
			 * 
			 * @param complex (Pointer to) a complex.
			 * @param dimension Subcomplex dimension.
			 * @param p _(Optional, default `0.5`)_ Bernoulli trial density.
			 * @param DEBUG _(Optional, default `false`)_ Are we debugging this Model?
			 */
			Bernoulli(
				ATEAMS::complexes::Complex<Z2>* complex,
				int dimension,
				float p=0.5,
				bool DEBUG=true
			);

			/**
			 * @brief Directly samples Bernoulli percolation at the desired dimension.
			 * 
			 * @param t Time step.
			 * @param state Tracks the state of the model. See @ref ATEAMS::statistics::Chain.
			 * @param options Multithreaded computing environment options.
			 * 
			 * @returns @ref ModelState with @ref ModelState::essential,
			 * @ref ModelState::includes, and @ref ModelState::t modified.
			 */
			ModelState<Z2,DenseVector> sample(
				int t, 
				ModelState<Z2,DenseVector>& state,
				ATEAMS::arithmetic::ThreadOptions& options
			) override;

			/** @brief Initialization; superfluous. */
			ModelState<Z2,DenseVector> initialize(ModelState<Z2,DenseVector>& state) override { return state; };

			/** @brief Initialization; superfluous. */
			ModelState<Z2,DenseVector> initialize(std::vector<Z2> c, ModelState<Z2,DenseVector>& state) override { return state; };

		private:
			std::mt19937 RNG;
			std::uniform_real_distribution<double> unituniform;
			std::vector<int> filtration;
	};
}

#include "ATEAMS++/models/Bernoulli.tpp"

#endif
