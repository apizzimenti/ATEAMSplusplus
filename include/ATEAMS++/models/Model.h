
#ifndef ATEAMS_MODELS_MODEL_H
#define ATEAMS_MODELS_MODEL_H

#include "ATEAMS++/complexes/Complex.h"

namespace ATEAMS {
	struct ModelParameters {
		int dimension;
	};


	class Model {
		public:
			Complex* complex;
			ModelParameters parameters;

			ZpVector sample();
	};
}

#endif
