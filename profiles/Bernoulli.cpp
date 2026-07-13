
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

using Model = models::Bernoulli;
using Parameters = models::ModelParameters;

using Structure = complexes::Cubical<Model::RingType>;
using State = models::Bernoulli::State;
using Chain = statistics::Chain<models::Bernoulli>;


int main(int argc, char* argv[]) {
	// cmd
	int SCALE = std::stoi(argv[1]);
	int TOPDIMENSION = std::stoi(argv[2]);
	int PLAQUETTEDIMENSION = std::stoi(argv[3]);
	int FIELD = std::stoi(argv[4]);
	int ITERATIONS = std::stoi(argv[5]);

	// Construct a cubical complex.
	vector<int> corners(TOPDIMENSION, SCALE);
	Structure COMPLEX(corners, true);

	// Parametrize + initialize the model.
	Parameters PARAMETERS;
	PARAMETERS.p = 0.5;
	PARAMETERS.dimension = PLAQUETTEDIMENSION;

	Model MODEL(&COMPLEX, PARAMETERS);
	Chain M(&MODEL, ITERATIONS);

	for (State STATE : M.simulate()) { }

	return 0;
}

