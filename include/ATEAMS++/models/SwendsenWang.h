
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
	 * @var SwendsenWang::parameters
	 * Model parameters.
	 * 
	 * @var SwendsenWang::state
	 * Model state.
	 * 
	 * @var SwendsenWang::field
	 * Finite field over which computations will be performened.
	 * 
	 * @var SwendsenWang::RNG
	 * Mersenne twister random number generator.
	 * 
	 * @var SwendsenWang::unituniform
	 * Uniform distribution over \f$(0,1)\f$.
	 * 
	 * @var SwendsenWang::intuniform
	 * Uniform distribution over the integers \f$[0,p)\f$, where \f$p\f$ is
	 * the order of the model's finite field.
	 * 
	 * @var SwendsenWang::kind
	 * Model name.
	 */
	template <typename T>
	class SwendsenWang: public Model<T,SparseVector> {
		public:
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
			 * @param options Multithreaded computing environment options.
			 * 
			 * @return The sample \f$f_{t+1}\f$.
			 */
			ModelState<SparseVector<T>> sample(int t, ATEAMS::arithmetic::ThreadOptions& options) override;

			/**
			 * @brief Initializes \f$f_0\f$ to uniform random elements of \f$\Z/p\Z\f$.
			 */
			void initialize() override;

			/**
			 * @brief Initializes \f$f_0 = c\f$.
			 */
			void initialize(SparseVector<T> c) override;

			ModelParameters parameters;
			ModelState<SparseVector<T>> state;
			const Field field;

			std::string kind = "SwendsenWang";

		private:
			std::mt19937 RNG;
			std::uniform_real_distribution<double> unituniform;
			std::uniform_int_distribution<int> intuniform;
			std::set<int> include;
	};

	#include "models/SwendsenWang.tpp"
}

#endif

