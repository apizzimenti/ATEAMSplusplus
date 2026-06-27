
#include <ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

// For ease-of-use.
typedef statistics::Chain<models::Bernoulli> Chain;

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
	models::BernoulliParameters params;
	params.p = 0.5;
	params.dimension = PLAQUETTEDIMENSION;

	models::Bernoulli G(&C, params);

	// Create the chain and data storage buckets.
	Chain M(&G, ITERATIONS);

	for (models::BernoulliState* state : M.simulate<models::BernoulliState>()) { }

	return 0;
}


