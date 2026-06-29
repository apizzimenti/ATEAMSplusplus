
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Parameters = models::BernoulliParameters;
using Model = models::Bernoulli;
using State = models::BernoulliState;
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
		PARAMETERS.p = 0.5;
		PARAMETERS.dimension = dimension/2;

		Model MODEL(&CUBICAL, PARAMETERS);
		Chain CHAIN(&MODEL, 10);

		for (State* state : CHAIN.simulate<State>()) {
			if (state->rank > ranks[dimension]) RESULT = FAIL;
		}
	}

	return RESULT;
}

