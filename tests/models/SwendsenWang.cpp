
#include <ATEAMS++/ATEAMS++.h>
#include <ATEAMS++/models/SwendsenWang.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <ATEAMS++/statistics/schedules.h>
#include <ATEAMS++/arithmetic/options.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {10,10,10,10};
	complexes::Cubical C(corners, true);

	models::SwendsenWangParameters params;
	params.field = 3;
	params.temperatureFunction = statistics::selfdual(params.field);
	params.dimension = 2;
	params.DEBUG = true;

	models::SwendsenWang SW(&C, params);
	SW.initialize();

	ATEAMS::arithmetic::ThreadOptions options;
	std::thread listener = options.spinUp();
	
	for (int i=0; i < 100; i++) {
		SW.sample(i, options);
	}

	options.spinDown(&listener);
	
	return 0;
}

