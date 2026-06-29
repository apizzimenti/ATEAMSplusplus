
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Complex = complexes::Cubical;
using Parameters = models::BernoulliParameters;
using Model = models::Bernoulli;
using State = models::BernoulliState;
using Chain = statistics::Chain<Model>;

int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	for (int dimension : DIMENSIONS) {
		vector<int> corners(dimension, 3);
		Complex COMPLEX(corners);

		Parameters PARAMETERS;
		PARAMETERS.p = 0.5;
		PARAMETERS.dimension = dimension/2;

		Model MODEL(&COMPLEX, PARAMETERS);
		Chain CHAIN(&MODEL, ITERATIONS);

		for (State* STATE : CHAIN.simulate<State>()) {
			if (STATE->rank > HOMOLOGICALRANK[dimension]) RESULT = FAIL;
		}
	}

	return RESULT;
}

