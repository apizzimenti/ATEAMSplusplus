
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Parameters = models::InvasionParameters;
using Model = models::Invasion;
using State = models::InvasionState;
using Chain = statistics::Chain<Model>;

int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	vector<int> dimensions{2,3,4};
	map<int,int> ranks{
		{2, 2},
		{3, 3},
		{4, 6}
	};

	for (int dimension : dimensions) {
		vector<int> corners(dimension, 3);
		complexes::Cubical CUBICAL(corners);

		Parameters PARAMETERS;
		PARAMETERS.field = FIELD;
		PARAMETERS.dimension = dimension/2;
		PARAMETERS.stoppingFunction = statistics::stopInvadingAt({1});

		Model MODEL(&CUBICAL, PARAMETERS);
		Chain CHAIN(&MODEL, 10);

		for (State* state : CHAIN.simulate<State>()) {
			if (state->essential.size() != ranks[dimension]) RESULT = FAIL;
		}
	}

	return RESULT;
}

