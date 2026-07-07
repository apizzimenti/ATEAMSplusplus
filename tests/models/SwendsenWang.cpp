
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Structure = complexes::Cubical<FINITE>;
using Model = models::SwendsenWang<FINITE>;
using State = models::ModelState<FINITE,SparseVector>;
using Chain = statistics::Chain<FINITE,SparseVector>;

int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	for (int dimension : DIMENSIONS) {
		vector<int> corners(dimension, 3);
		Structure CUBICAL(corners);

		models::ModelParameters PARAMETERS;
		PARAMETERS.field = FIELD;
		PARAMETERS.dimension = dimension/2;
		PARAMETERS.temperatureFunction = statistics::selfdual(PARAMETERS.field);

		Model MODEL(&CUBICAL, PARAMETERS);
		Chain CHAIN(&MODEL, ITERATIONS);

		for (State state : CHAIN.simulate()) {
			// Figure out which cells were excluded; on these cells, the cochain
			// can evaluate to anything. We just want the ones that evaluate to
			// 0.
			vector<int> unsatisfied = statistics::unsatisfied(
				&CUBICAL, state.cochain, MODEL.field, PARAMETERS.dimension
			);

			ZpMatrix REDUCED = MODEL.complex->Coboundary.Matrices[PARAMETERS.dimension];
			for (auto u : unsatisfied) REDUCED[u].zero();
			REDUCED.clear_zero_row();
			REDUCED.compress();

			if (!inKernel(REDUCED, state.cochain, MODEL.field)) {
				RESULT = FAIL;
			}
		}
	}

	return RESULT;
}

