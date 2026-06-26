
#ifndef ATEAMS_MODELS_BERNOULLI_H
#define ATEAMS_MODELS_BERNOULLI_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/persistence.h"

#include <random>
#include <set>

namespace ATEAMS::models {
	/**
	 * @struct BernoulliParameters
	 * @brief Parameters for Bernoulli percolation.
	 * 
	 * @var BernoulliParameters::p
	 * Bernoulli trial density.
	 */
	struct BernoulliParameters : ModelParameters {
		double p = 0.5;
	};

	/**
	 * @struct BernoulliState
	 * @brief Stores the current state of Bernoulli percolation.
	 * 
	 * @var BernoulliState::includes
	 * Indices of included \f$d\f$-cells.
	 * 
	 * @var BernoulliState::essential
	 * Times at which homological percolation occurred (i.e. giant cycles appeared).
	 * 
	 * @var BernoulliState::rank
	 * Rank of the \f$d\f$th persistent homology group \f$PH_d(P)\f$.
	 */
	struct BernoulliState : ModelState {
		std::vector<int> includes;
		std::vector<int> essential;
		int rank;
	};

	/**
	 * @class Bernoulli
	 * @brief Implements Bernoulli percolation.
	 * 
	 * @var Bernoulli::parameters
	 * Parameters for Bernoulli percolation; includes standard ModelParameters, and
	 * 	Bernoulli trial density \f$0 \leq p \leq 1\f$.
	 * 
	 * @var Bernoulli::state
	 * Records state.
	 * 
	 * @var Bernoulli::kind
	 * Model name.
	 */
	class Bernoulli : Model<std::vector<int>> {
		public:
			/**
			 * @brief Constructor.
			 * @param complex (Pointer to) a complex.
			 * @param parameters Model parameters.
			 */
			Bernoulli(ATEAMS::complexes::Complex* complex, BernoulliParameters parameters);

			/**
			 * @brief Directly samples Bernoulli percolation at the desired dimension.
			 * @param t Time step.
			 * @param options Multithreaded computing environment options.
			 * 
			 * @return List of indices of included \f$d\f$-cells.
			 */
			std::vector<int> sample(int t, ATEAMS::arithmetic::ThreadOptions& options) override;

			/** @brief Initialization; superfluous. */
			void initialize() override { };

			/** @brief Initialization; superfluous. */
			void initialize(std::vector<int> c) override { };

			BernoulliParameters parameters;
			BernoulliState state;

			std::string kind = "Bernoulli";

		private:
			std::mt19937 RNG;
			std::uniform_real_distribution<double> unituniform;
			const Zp field;
			
			std::vector<int> filtration;
	};
}

#endif
