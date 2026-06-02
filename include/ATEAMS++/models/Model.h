
#ifndef ATEAMS_MODELS_MODEL_H
#define ATEAMS_MODELS_MODEL_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/complexes/Complex.h"

namespace ATEAMS {
	struct ModelParameters {
		int dimension;
	};

	struct ModelState { };


	class Model {
		public:
			Complex* complex;
			ModelParameters parameters;
			ModelState state;

			virtual ZpVector sample(int t) = 0;
			virtual void initial() = 0;
			virtual void initial(ZpVector c) = 0;
	};
}

#endif
