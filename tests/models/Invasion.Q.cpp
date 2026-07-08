
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Structure = complexes::Cubical<FINITE>;
using Model = models::Invasion;
using State = models::ModelState<FINITE,DenseVector>;
using Chain = statistics::Chain<FINITE,DenseVector>;

int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	for (int dimension : DIMENSIONS) {
		vector<int> corners(dimension, 3);
		Structure COMPLEX(corners);

		models::ModelParameters PARAMETERS;
		PARAMETERS.field = FIELD;
		PARAMETERS.dimension = dimension/2;
		PARAMETERS.stoppingFunction = statistics::stopInvadingAt({1});

		Model MODEL(&COMPLEX, PARAMETERS);
		Chain CHAIN(&MODEL, ITERATIONS);

		for (State state : CHAIN.simulate()) {
			if (state.essential.size() != HOMOLOGICALRANK[dimension]) RESULT = FAIL;
		}
	}

	return RESULT;
}

