
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
	 */
	template <typename T>
	class SwendsenWang: public Model<T,SparseVector> {
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
			SwendsenWang(ATEAMS::complexes::Complex<T>* complex, ModelParameters parameters);

			/**
			 * @brief Alternates between sampling the conditionals \f$P_t \sim \PK(- \mid f_t)\f$
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
			ModelState<T,SparseVector> sample(int t, ModelState<T,SparseVector>& state, ATEAMS::arithmetic::ThreadOptions& options) override;

			/**
			 * @brief Initializes \f$f_0\f$ to uniform random element of \f$\Z/p\Z\f$
			 * or \f$ \Q \f$.
			 * 
			 * @param state State.
			 * 
			 * @returns Modified State.
			 */
			ModelState<T,SparseVector> initialize(ModelState<T,SparseVector>& state) override;

			/**
			 * @brief Initializes \f$f_0 = c\f$.
			 * 
			 * @param c Cochain to which \f$f_0\f$ is initialized.
			 * @param state State.
			 * 
			 * @returns Modified State.
			 */
			ModelState<T,SparseVector> initialize(SparseVector<T> c, ModelState<T,SparseVector>& state) override;

		private:
			std::mt19937 RNG;
			std::uniform_real_distribution<double> unituniform;
			std::uniform_int_distribution<int> intuniform;
			std::set<int> include;
	};
}

#include "ATEAMS++/models/SwendsenWang.tpp"

#endif

