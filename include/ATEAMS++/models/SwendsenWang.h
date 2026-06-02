
#ifndef ATEAMS_MODELS_SWENDSENWANG_H
#define ATEAMS_MODELS_SWENDSENWANG_H

#include "ATEAMS++/models/Model.h"

namespace ATEAMS {
	struct SwendsenWangParameters : ModelParameters {
		int field;
		double temperatureFunction();
		Zp F;
	};


	class SwendsenWang: public Model {
		public:
			SwendsenWang(Complex* complex, SwendsenWangParameters parameters);
			ZpVector sample();

			SwendsenWangParameters parameters;

		private:
			void delegateComputation();
	};
}

#endif

