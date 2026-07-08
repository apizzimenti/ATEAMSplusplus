
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
		 * @class Glauber
		 * @brief Implements the Glauber dynamics algorithm.
		 * 
		 * @var Glauber::parameters
		 * Model parameters.
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
		template <typename T=ATEAMS::ff>
		class Glauber: public Model<T,SparseVector> {
			public:
				/**
				 * @brief Constructor.
				 * 
				 * @param complex (Pointer to) a complex.
				 * @param parameters Model parameters.
				 */
				Glauber(ATEAMS::complexes::Complex<T>* complex, ModelParameters parameters);

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
				ModelState<T,SparseVector> sample(int t, ModelState<T,SparseVector>& state, ATEAMS::arithmetic::ThreadOptions& options) override;

				/**
				 * @brief Initializes \f$f_0\f$ to uniform random element of \f$\Z/p\Z\f$
				 * or \f$ \Q \f$.
				 * 
				 * @param state Model state.
				 * 
				 * @returns Model state.
				 */
				ModelState<T,SparseVector> initialize(ModelState<T,SparseVector>& state) override;

				/**
				 * @brief Initializes \f$f_0 = c\f$.
				 * 
				 * @param c Cochain to which \f$f_0\f$ is initialized.
				 * @param state Model state.
				 * 
				 * @returns Model state.
				 */
				ModelState<T,SparseVector> initialize(SparseVector<T> c, ModelState<T,SparseVector>& state) override;
				
				ModelParameters parameters;
				const Field field;

				std::string kind = "Glauber";

			private:
				std::mt19937 RNG;
				std::uniform_real_distribution<double> unituniform;
				std::uniform_int_distribution<int> intuniform;
				std::uniform_int_distribution<int> indexuniform;
		};
	}
}

#include "ATEAMS++/models/Glauber.tpp"

#endif

