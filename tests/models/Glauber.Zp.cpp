
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Model = models::Glauber<Zp>;
using Parameters = models::ModelParameters;

using Structure = complexes::Cubical<Model::RingType>;
using State = Model::State;
using Chain = statistics::Chain<Model>;

int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	for (int dimension : DIMENSIONS) {
		vector<int> corners(dimension, 3);
		Structure CUBICAL(corners);

		Model::RingType RR(FIELD);

		Parameters PARAMETERS;
		PARAMETERS.coefficients = &RR;
		PARAMETERS.dimension = dimension/2;
		PARAMETERS.temperatureFunction = statistics::selfdual(PARAMETERS.coefficients);
		PARAMETERS.DEBUG = true;

		Model MODEL(&CUBICAL, PARAMETERS);
		Chain CHAIN(&MODEL, ITERATIONS);

		for (State state : CHAIN.simulate()) {
			// Figure out which cells were excluded; on these cells, the cochain
			// can evaluate to anything. We just want the ones that evaluate to
			// 0.
			vector<int> unsatisfied = statistics::unsatisfied<Model::RingType>(
				&CUBICAL, state.cochain, MODEL.coefficients, PARAMETERS.dimension
			);

			SparseMatrix<Model::RingType> REDUCED = MODEL.complex->Coboundary.Matrices[PARAMETERS.dimension];
			for (auto u : unsatisfied) REDUCED[u].zero();
			REDUCED.clear_zero_row();
			REDUCED.compress();

			if (!inKernel<Model::RingType>(REDUCED, state.cochain, MODEL.coefficients, PARAMETERS.DEBUG)) {
				RESULT = FAIL;
			}
		}
	}

	return RESULT;
}

