
#include <ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

// For ease-of-use.
typedef statistics::Chain<models::InvadedCluster> Chain;

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
	complexes::Cubical C(corners, true);

	// Parametrize + initialize the model.
	models::InvadedClusterParameters params;
	params.dimension = PLAQUETTEDIMENSION;
	params.field = FIELD;
	params.stoppingFunction = statistics::stopInvadingAt(stopat[TOPDIMENSION]);

	models::InvadedCluster G(&C, params);

	// Create the chain and data storage buckets.
	Chain M(&G, ITERATIONS);

	for (models::InvadedClusterState* state : M.simulate<models::InvadedClusterState>()) { }

	return 0;
}


