
#include <ATEAMS++/models/SwendsenWang.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <ATEAMS++/statistics/schedules.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {10,10,10,10};
	Cubical C(corners, true);

	SwendsenWangParameters params;
	params.field = 3;
	params.temperatureFunction = ATEAMS::critical(params.field);
	params.dimension = 2;
	params._DEBUG = true;

	SwendsenWang SW(&C, params);
	SW.initial();

	for (int i=0; i < 200; i++) {
		SW.sample(i);
	}
	
	return 0;
}

