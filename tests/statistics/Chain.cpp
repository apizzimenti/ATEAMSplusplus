
#include <ATEAMS++/ATEAMS++.h>
#include <ATEAMS++/statistics/Chain.h>
#include <ATEAMS++/models/InvadedCluster.h>
#include <ATEAMS++/complexes/Cubical.h>

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

	using Chain = statistics::Chain<models::InvadedCluster>;
	Chain M(&SW, 100);

	for (models::InvadedClusterState* state : M.simulate<models::InvadedClusterState>()) {
		// Typecasting.
		cout << state->cochain.size() << endl;
	}
	
	return 0;
}


