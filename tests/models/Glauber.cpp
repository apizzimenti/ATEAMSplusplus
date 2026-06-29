
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Complex = complexes::Cubical;
using Parameters = models::GlauberParameters;
using Model = models::Glauber;
using State = models::GlauberState;
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
		PARAMETERS.temperatureFunction = statistics::constant(PARAMETERS.field);

		Model MODEL(&COMPLEX, PARAMETERS);
		Chain CHAIN(&MODEL, ITERATIONS);

		for (State* state : CHAIN.simulate<State>()) {
			// Figure out which cells were excluded; on these cells, the cochain
			// can evaluate to anything. We just want the ones that evaluate to
			// 0.
			vector<int> unsatisfied = statistics::unsatisfied(
				&COMPLEX, state->cochain, MODEL.field, PARAMETERS.dimension
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

