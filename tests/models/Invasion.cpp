
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Complex = complexes::Cubical;
using Parameters = models::InvasionParameters;
using Model = models::Invasion;
using State = models::InvasionState;
using Chain = statistics::Chain<Model>;

int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	for (int dimension : DIMENSIONS) {
		vector<int> corners(dimension, 3);
		Complex COMPLEX(corners);

		Parameters PARAMETERS;
		PARAMETERS.field = FIELD;
		PARAMETERS.dimension = dimension/2;
		PARAMETERS.stoppingFunction = statistics::stopInvadingAt({1});

		Model MODEL(&COMPLEX, PARAMETERS);
		Chain CHAIN(&MODEL, ITERATIONS);

		for (State* state : CHAIN.simulate<State>()) {
			if (state->essential.size() != HOMOLOGICALRANK[dimension]) RESULT = FAIL;
		}
	}

	return RESULT;
}

