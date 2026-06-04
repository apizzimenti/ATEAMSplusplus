
#include <ATEAMS++/util.h>
#include <ATEAMS++/common.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <iostream>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {3,3,3,3};
	complexes::Cubical C(corners, true);
	ATEAMS::Zp F(SparseRREF::FIELD_Fp, 2);

	C.constructBoundaryMatrices(F);

	for (int d=0; d < corners.size()+1; d++) {
		printvector(C.Boundary.breaks[d]);
	}
	
	return 0;
}
