
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {6,6,6,6};
	complexes::Cubical C(corners, true);

	models::GlauberParameters params;
	params.field = 3;
	params.temperatureFunction = statistics::selfdual(params.field);
	params.dimension = 2;
	params.DEBUG = false;

	models::Glauber Glauber(&C, params);

	using Chain = statistics::Chain<models::Glauber>;
	Chain M(&Glauber, 1000);

	for (models::ModelState* state : M.simulate()) { }

	return 0;
}

