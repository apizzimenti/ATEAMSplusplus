
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {6,6,6,6};
	complexes::Cubical C(corners, true);

	models::InvasionPercolationParameters params;
	params.stoppingFunction = arithmetic::stopInvadingAt({3,4});
	params.dimension = 2;

	models::InvasionPercolation percolation(&C, params);

	using Chain = statistics::Chain<models::InvasionPercolation>;
	Chain M(&percolation, 1000);

	for (models::ModelState* state : M.simulate()) { };
	
	return 0;
}

