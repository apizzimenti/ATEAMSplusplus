
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
		 */
		template <typename T=ATEAMS::ff>
		class Glauber: public Model<T,SparseVector> {
			public:
				/** Coefficients are of type @ref ATEAMS::ff or @ref ATEAMS::rational. */
				typedef T dt;

				/** Vectors are of type @ref ATEAMS::SparseVector. */
				template <typename R>
				using st = SparseVector<R>;
				
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

