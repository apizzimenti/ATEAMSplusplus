
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
	 * @var Bernoulli::parameters
	 * Parameters for Bernoulli percolation; includes standard ModelParameters, and
	 * 	Bernoulli trial density \f$0 \leq p \leq 1\f$.
	 * 
	 * @var Bernoulli::kind
	 * Model name.
	 */
	template <typename T>
	class Bernoulli : Model<T,std::vector> {
		public:
			ModelParameters parameters;
			std::string kind = "Bernoulli"

			/**
			 * @brief Constructor.
			 * @param complex (Pointer to) a complex.
			 * @param parameters Model parameters.
			 */
			Bernoulli(ATEAMS::complexes::Complex<T>* complex, BernoulliParameters parameters);

			/**
			 * @brief Directly samples Bernoulli percolation at the desired dimension.
			 * @param t Time step.
			 * @param options Multithreaded computing environment options.
			 * 
			 * @return List of indices of included \f$d\f$-cells.
			 */
			ModelState<T,std::vector> sample(int t, ModelState<T,std::vector>& state, ATEAMS::arithmetic::ThreadOptions& options) override;

			/** @brief Initialization; superfluous. */
			ModelState<T,std::vector> initialize(ModelState<T,std::vector>& state) override { };

			/** @brief Initialization; superfluous. */
			ModelState<T,std::vector> initialize(std::vector<int> c, ModelState<T,std::vector>& state) override { };

		private:
			std::mt19937 RNG;
			std::uniform_real_distribution<double> unituniform;
			const Field field;
			
			std::vector<int> filtration;
	};
}

#endif
