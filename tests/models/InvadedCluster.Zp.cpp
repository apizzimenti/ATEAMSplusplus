
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

using Model = models::InvadedCluster<FINITE>;
using Parameters = models::ModelParameters;

using Structure = complexes::Cubical<Model::dt>;
using State = models::ModelState<Model::dt,Model::st>;
using Chain = statistics::Chain<Model::dt,Model::st>;

int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	for (int dimension : DIMENSIONS) {
		vector<int> corners(dimension, 3);
		Structure COMPLEX(corners);

		Parameters PARAMETERS;
		PARAMETERS.field = FIELD;
		PARAMETERS.dimension = dimension/2;
		PARAMETERS.stoppingFunction = statistics::stopInvadingAt({dimension/2});
		PARAMETERS.DEBUG = true;

		Model MODEL(&COMPLEX, PARAMETERS);
		Chain CHAIN(&MODEL, ITERATIONS);

		for (State state : CHAIN.simulate()) {
			// Figure out which cells were excluded; on these cells, the cochain
			// can evaluate to anything. We just want the ones that evaluate to
			// 0.
			vector<int> unsatisfied = statistics::unsatisfied<Model::dt>(
				&COMPLEX, state.cochain, MODEL.field, PARAMETERS.dimension
			);

			SparseMatrix<Model::dt> REDUCED = MODEL.complex->Coboundary.Matrices[PARAMETERS.dimension];
			for (auto u : unsatisfied) REDUCED[u].zero();
			REDUCED.clear_zero_row();
			REDUCED.compress();

			if (!inKernel<Model::dt>(REDUCED, state.cochain, MODEL.field, PARAMETERS.DEBUG)) {
				RESULT = FAIL;
			}
		}
	}

	return RESULT;
}

