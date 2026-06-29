
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Parameters = models::InvadedClusterParameters;
using Model = models::InvadedCluster;
using State = models::InvadedClusterState;
using Chain = statistics::Chain<Model>;

int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	vector<int> dimensions{2,3,4};

	for (int dimension : dimensions) {
		vector<int> corners(dimension, 3);
		complexes::Cubical CUBICAL(corners);

		Parameters PARAMETERS;
		PARAMETERS.field = FIELD;
		PARAMETERS.dimension = dimension/2;
		PARAMETERS.stoppingFunction = statistics::stopInvadingAt({dimension/2});

		Model MODEL(&CUBICAL, PARAMETERS);
		Chain CHAIN(&MODEL, 10);

		for (State* state : CHAIN.simulate<State>()) {
			// Figure out which cells were excluded; on these cells, the cochain
			// can evaluate to anything. We just want the ones that evaluate to
			// 0.
			vector<int> unsatisfied = statistics::unsatisfied(
				&CUBICAL, state->cochain, MODEL.field, PARAMETERS.dimension
			);

			ZpMatrix REDUCED = MODEL.complex->Coboundary.Matrices[PARAMETERS.dimension];
			for (auto u : unsatisfied) REDUCED[u].zero();
			REDUCED.clear_zero_row();
			REDUCED.compress();

			if (!inKernel(REDUCED, state->cochain, MODEL.field)) {
				RESULT = FAIL;
			}
		}
	}

	return RESULT;
}

