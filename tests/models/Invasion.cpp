
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {6,6,6,6};
	complexes::Cubical C(corners, true);

	models::InvasionParameters params;
	params.stoppingFunction = arithmetic::stopInvadingAt({3,4});
	params.dimension = 2;

	models::Invasion percolation(&C, params);

	using Chain = statistics::Chain<models::Invasion>;
	Chain M(&percolation, 1000);

	for (models::ModelState* state : M.simulate()) { };
	
	return 0;
}

