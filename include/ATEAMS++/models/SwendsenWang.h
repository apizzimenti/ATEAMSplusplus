
#ifndef ATEAMS_MODELS_SWENDSENWANG_H
#define ATEAMS_MODELS_SWENDSENWANG_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Model.h"

#include <random>

namespace ATEAMS {
	struct SwendsenWangParameters : ModelParameters {
		int field = 2;
		std::function<double(int)> temperatureFunction;
	};

	struct SwendsenWangState : ModelState {
		ZpVector cochain;
		ZpVector includes;
	};


	class SwendsenWang: public Model {
		public:
			SwendsenWang(Complex* complex, SwendsenWangParameters parameters);

			ZpVector sample(int t) override;

			void initial() override;
			void initial(ZpVector c) override;

			SwendsenWangParameters parameters;
			SwendsenWangState state;
			const Zp field;

		private:
			void delegateComputation();

			// random number generation
			std::mt19937 RNG;
			std::uniform_real_distribution<double> unituniform;
			std::uniform_int_distribution<int> intuniform;
	};
}

#endif

