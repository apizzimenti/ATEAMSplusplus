
#include <ATEAMS++/models/SwendsenWang.h>
#include <ATEAMS++/complexes/Cubical.h>

using namespace ATEAMS;
using namespace std;

double temp(int t) {
	return (double)0;
}

int main() {
	vector<int> corners = {3,3};
	Cubical C(corners, true);

	SwendsenWangParameters params;
	params.temperatureFunction = temp;
	params.field = 3;
	params.dimension = 1;

	SwendsenWang SW(&C, params);
	SW.initial();
	SW.sample(0);
	
	return 0;
}

