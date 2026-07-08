
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Structure = complexes::Cubical<FINITE>;
using Model = models::Bernoulli;
using State = models::ModelState<FINITE,DenseVector>;
using Chain = statistics::Chain<FINITE,DenseVector>;

int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	for (int dimension : DIMENSIONS) {
		vector<int> corners(dimension, 3);
		Structure COMPLEX(corners);

		models::ModelParameters PARAMETERS;
		PARAMETERS.p = 0.5;
		PARAMETERS.dimension = dimension/2;

		Model MODEL(&COMPLEX, PARAMETERS);
		Chain CHAIN(&MODEL, ITERATIONS);

		for (State STATE : CHAIN.simulate()) {
			if (STATE.rank > HOMOLOGICALRANK[dimension]) RESULT = FAIL;
		}
	}

	return RESULT;
}

