
#include <ATEAMS++/util.h>
#include <ATEAMS++/common.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <ATEAMS++/arithmetic/persistence.h>
#include <iostream>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {3,3};
	complexes::Cubical C(corners, true);
	ATEAMS::Zp F(SparseRREF::FIELD_Fp, 2);

	C.constructFlatBoundaryMatrix();
	C.constructFullBoundaryMatrix(F);

	vector<int> filtration(C.size());
	std::iota(begin(filtration), end(filtration), 0);

	filtration[16] = 17;
	filtration[17] = 16;

	C.constructFlatBoundaryMatrix();
	C.constructFullBoundaryMatrix(F);

	vector<int> essential = arithmetic::PHATPersistence(&C, filtration, 1);
	vector<int> essentialmin = arithmetic::TwistPersistence(&C, filtration, F, 1);

	printvector<int>(essential);
	printvector<int>(essentialmin);
	
	return 0;
}
