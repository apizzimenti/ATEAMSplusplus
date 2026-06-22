
#ifndef ATEAMS_MODELS_GLAUBER_H
#define ATEAMS_MODELS_GLAUBER_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"
#include "ATEAMS++/arithmetic/options.h"

#include <random>
#include <string>

namespace ATEAMS {
	namespace models {
		/**
		 * @struct GlauberParameters
		 * @brief Parameters relevant to simulating the Glauber dynamics Markov chain.
		 * 
		 * @var GlauberParameters::field
		 * Order \f$p\f$ for the finite field \f$\Z / p \Z\f$.
		 * 
		 * @var GlauberParameters::dimension
		 * Complex dimension for which we're sampling.
		 * 
		 * @var GlauberParameters::temperatureFunction
		 * Function that consumes a time-step (integer) parameter \f$t\f$, and returns
		 * a temperature parameter \f$\beta \in [0, \infty)\f$ that is used to compute
		 * the probability \f$p = 1-e^{-\beta}\f$ of including any individual
		 * \f$(d-1)\f$-cell in the sampled subcomplex.
		 */
		struct GlauberParameters : ModelParameters {
			int field = 2;
			int dimension;
			std::function<double(int)> temperatureFunction;
		};

		/**
		 * @struct GlauberState
		 * @brief Stores the current state of the Glauber dynamics Markov chain.
		 * 
		 * @var GlauberState::cochain
		 * A (sparse) vector representing the cochain \f$f_t \in C^{d-1}(X; \Z/p\Z)\f$
		 * (i.e. a linear map from \f$(d-1)\f$-cells to \f$\Z/p\Z\f$)
		 * sampled at time \f$t\f$.
		 * 
		 * @var GlauberState::energy
		 * An integer representing the total energy of the cochain (i.e. the Hamiltonian
		 * \f$\H(f_t)\f$) indicating the (negative) number of \f$d\f$-cells \f$x\f$ on which
		 * \f$(\delta^{d-1}f)(x) > 0\f$.
		 */
		struct GlauberState : ModelState {
			ZpVector cochain;
			int energy;
		};

		/**
		 * @class Glauber
		 * @brief Implements the Glauber dynamics algorithm.
		 * 
		 * @var Glauber::parameters
		 * Model parameters.
		 * 
		 * @var Glauber::state
		 * Model state.
		 * 
		 * @var Glauber::field
		 * Finite field over which computations will be performened.
		 * 
		 * @var Glauber::RNG
		 * Mersenne twister random number generator.
		 * 
		 * @var Glauber::unituniform
		 * Uniform distribution over \f$(0,1)\f$.
		 * 
		 * @var Glauber::intuniform
		 * Uniform distribution over the integers \f$[0,p)\f$, where \f$p\f$ is
		 * the order of the model's finite field.
		 * 
		 * @var Glauber::indexuniform
		 * Uniform distribution over the integers \f$[0,N)\f$, where \f$N\f$ is
		 * the number of \f$(d-1)\f$-cells in the complex.
		 * 
		 * @var Glauber::kind
		 * Model name.
		 */
		class Glauber: public Model<ZpVector> {
			public:
				/**
				 * @brief Constructor.
				 * 
				 * @param complex (Pointer to) a complex.
				 * @param parameters Model parameters.
				 */
				Glauber(ATEAMS::complexes::Complex* complex, GlauberParameters parameters);

				/**
				 * @brief Switches the spin of a uniform random cell.
				 * 
				 * @param t Time step.
				 * @param options Multithreaded computing environment options.
				 * 
				 * @return The sample \f$f_{t+1}\f$.
				 */
				ZpVector sample(int t, ATEAMS::arithmetic::ThreadOptions& options) override;

				/**
				 * @brief Initializes \f$f_0\f$ to uniform random elements of \f$\Z/p\Z\f$.
				 */
				void initialize() override;

				/**
				 * @brief Initializes \f$f_0 = c\f$.
				 */
				void initialize(ZpVector c) override;

				GlauberParameters parameters;
				GlauberState state;
				const Zp field;

				std::string kind = "Glauber";

			private:
				std::mt19937 RNG;
				std::uniform_real_distribution<double> unituniform;
				std::uniform_int_distribution<int> intuniform;
				std::uniform_int_distribution<int> indexuniform;
		};
	}
}

#endif

