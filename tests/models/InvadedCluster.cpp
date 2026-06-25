
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {6,6,6,6};
	complexes::Cubical C(corners, true);

	models::InvadedClusterParameters params;
	params.field = 3;
	params.stoppingFunction = arithmetic::stopInvadingAt({3,4});
	params.dimension = 2;
	params.DEBUG = false;

	models::InvadedCluster invasion(&C, params);

	using Chain = statistics::Chain<models::InvadedCluster>;
	Chain M(&invasion, 1000);

	for (models::ModelState* state : M.simulate()) { };
	
	return 0;
}

