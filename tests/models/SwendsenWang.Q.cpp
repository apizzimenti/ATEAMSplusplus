
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Model = models::SwendsenWang<Q>;
using Parameters = models::ModelParameters;

using Structure = complexes::Cubical<Q>;
using State = models::ModelState<Q,SparseVector>;
using Chain = statistics::Chain<Q,SparseVector>;


int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	for (int dimension : DIMENSIONS) {
		vector<int> corners(dimension, 3);
		Structure CUBICAL(corners);

		Q QQ;

		Parameters PARAMETERS;
		PARAMETERS.coefficients = &QQ;
		PARAMETERS.dimension = dimension/2;
		PARAMETERS.temperatureFunction = statistics::constant(1/2);
		PARAMETERS.DEBUG = true;

		Model MODEL(&CUBICAL, PARAMETERS);
		Chain CHAIN(&MODEL, ITERATIONS);

		for (State state : CHAIN.simulate()) {
			// Figure out which cells were excluded; on these cells, the cochain
			// can evaluate to anything. We just want the ones that evaluate to
			// 0.
			vector<int> unsatisfied = statistics::unsatisfied<Q>(
				&CUBICAL, state.cochain, MODEL.coefficients, PARAMETERS.dimension
			);

			SparseMatrix<Q> REDUCED = MODEL.complex->Coboundary.Matrices[PARAMETERS.dimension];
			for (auto u : unsatisfied) REDUCED[u].zero();
			REDUCED.clear_zero_row();
			REDUCED.compress();

			if (!inKernel<Q>(REDUCED, state.cochain, MODEL.coefficients)) {
				RESULT = FAIL;
			}
		}
	}

	return RESULT;
}

