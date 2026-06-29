
#include <ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

using Complex = complexes::Cubical;
using Parameters = models::InvadedClusterParameters;
using Model = models::InvadedCluster;
using State = models::InvadedClusterState;
using Chain = statistics::Chain<Model>;

std::map<int,vector<int>> stopat = {
	{2,{1}},
	{4,{3,4}}
};

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
	params.stoppingFunction = statistics::stopInvadingAt(stopat[TOPDIMENSION]);

	Model MODEL(&COMPLEX, PARAMETERS);
	Chain M(&MODEL, ITERATIONS);

	for (State* STATE : M.simulate<State>()) { }
	
	return 0;
}


