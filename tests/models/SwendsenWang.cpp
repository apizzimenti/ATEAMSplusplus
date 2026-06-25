
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {6,6,6,6};
	complexes::Cubical C(corners, true);

	models::SwendsenWangParameters params;
	params.field = 3;
	params.temperatureFunction = statistics::selfdual(params.field);
	params.dimension = 2;

	models::SwendsenWang SW(&C, params);
	
	using Chain = statistics::Chain<models::SwendsenWang>;
	Chain M(&SW, 1000);

	for (models::ModelState* state : M.simulate()) { };
	
	return 0;
}

