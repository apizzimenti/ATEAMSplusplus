
#include <ATEAMS++/ATEAMS++.h>
#include <ATEAMS++/models/SwendsenWang.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <ATEAMS++/statistics/schedules.h>
#include <ATEAMS++/arithmetic/options.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {4,4,4,4};
	complexes::Cubical C(corners, true);

	models::SwendsenWangParameters params;
	params.field = 2;
	params.temperatureFunction = statistics::selfdual(params.field);
	params.dimension = 2;

	models::SwendsenWang SW(&C, params);
	SW.initialize();

	ATEAMS::arithmetic::ThreadOptions options;
	std::thread listener = options.spinUp();

	int row = SW.complex->Coboundary.Matrices[2].nrow;

	int k = 0;

	for (int t=0; t < row; t++) {
		k += SW.complex->Coboundary.Matrices[2][t].size();
	}

	options.spinDown(&listener);
	
	return 0;
}

