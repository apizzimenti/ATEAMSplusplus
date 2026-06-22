
#include <ATEAMS++/ATEAMS++.h>
#include <ATEAMS++/models/InvadedCluster.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <ATEAMS++/arithmetic/persistence.h>
#include <ATEAMS++/arithmetic/options.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {5,5,5,5};
	complexes::Cubical C(corners, true);

	models::InvadedClusterParameters params;
	params.field = 3;
	params.stoppingFunction = arithmetic::stopInvadingAt({3,4});
	params.dimension = 2;
	params.DEBUG = false;

	models::InvadedCluster SW(&C, params);
	SW.initialize();

	ATEAMS::arithmetic::ThreadOptions options;
	std::thread listener = options.spinUp();

	int N = 1000;

	for (int t=0; t < N; t++) {
		SW.sample(t, options);
	}

	options.spinDown(&listener);
	
	return 0;
}

