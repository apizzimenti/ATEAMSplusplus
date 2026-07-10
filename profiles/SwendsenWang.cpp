
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

using Model = models::SwendsenWang<Zp>;
using Parameters = models::ModelParameters;

using Structure = complexes::Cubical<Model::RingType>;
using State = models::ModelState<Model::RingType,Model::VectorType>;
using Chain = statistics::Chain<Model::RingType,Model::VectorType>;


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

	Model::RingType RR(FIELD);

	// Parametrize + initialize the model.
	Parameters PARAMETERS;
	PARAMETERS.dimension = PLAQUETTEDIMENSION;
	PARAMETERS.coefficients = &RR;
	PARAMETERS.temperatureFunction = statistics::selfdual(&RR);

	Model MODEL(&COMPLEX, PARAMETERS);
	Chain M(&MODEL, ITERATIONS);

	for (State STATE : M.simulate()) { }
	
	return 0;
}


