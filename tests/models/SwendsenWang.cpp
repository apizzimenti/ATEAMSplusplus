
#include <ATEAMS++/models/SwendsenWang.h>
#include <ATEAMS++/complexes/Cubical.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {3,3};
	Cubical C(corners, true);

	SwendsenWangParameters params;
	params.field = 3;

	SwendsenWang SW(&C, params);

	for (int d=0; d < SW.complex->Boundary.Matrices.size(); d++) {
		printSparseRREFmat(SW.complex->Boundary.Matrices[d]);
	}
	
	return 0;
}

