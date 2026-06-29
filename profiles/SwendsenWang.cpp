
#include <ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

using Complex = complexes::Cubical;
using Parameters = models::SwendsenWangParameters;
using Model = models::SwendsenWang;
using State = models::SwendsenWangState;
using Chain = statistics::Chain<Model>;

int main(int argc, char* argv[]) {
	// cmd
	int SCALE = std::stoi(argv[1]);
	int TOPDIMENSION = std::stoi(argv[2]);
	int PLAQUETTEDIMENSION = std::stoi(argv[3]);
	int FIELD = std::stoi(argv[4]);
	int ITERATIONS = std::stoi(argv[5]);

	// Construct a cubical complex.
	vector<int> corners(TOPDIMENSION, SCALE);
	Complex COMPLEX(corners, true);

	// Parametrize + initialize the model.
	Parameters PARAMETERS;
	PARAMETERS.dimension = PLAQUETTEDIMENSION;
	PARAMETERS.field = FIELD;
	PARAMETERS.temperatureFunction = statistics::selfdual(FIELD);

	Model MODEL(&COMPLEX, PARAMETERS);
	Chain M(&MODEL, ITERATIONS);

	for (State* STATE : M.simulate<State>()) { }
	return 0;
}


