
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {3,3,3,3};
	complexes::Cubical C(corners, true);

	models::BernoulliParameters params;
	params.p = 0.5;
	params.dimension = 2;

	models::Bernoulli percolation(&C, params);

	using Chain = statistics::Chain<models::Bernoulli>;
	Chain M(&percolation, 10);

	for (models::ModelState* state : M.simulate()) { }

	return 0;
}

