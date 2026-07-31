
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

#include <random>
#include <SparseRREF/sparse_mat.h>

using namespace ATEAMS;
using namespace std;


int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	// Random device, uniform distribution.
	std::random_device rd;
	std::mt19937 RNG(rd());
	std::uniform_int_distribution<int> intuniform(0, FIELD);

	// Construct a Complex (doesn't matter which).
	complexes::Cubical<Zp> C({3,3,3,3});
	Zp ZZ(FIELD);
	C.constructBoundaryMatrices(&ZZ);

	// Construct arithmetic resources.
	arithmetic::ComputeResources resources;
	std::thread listener = resources.spinUp();

	// Include all the rows and perform the computation.
	set<size_t> exclude;

	SparseVector<Zp> sample = arithmetic::submatrixKernelSample<Zp>(
		C.Coboundary.Matrices[2],
		&ZZ,
		exclude,
		intuniform,
		RNG,
		resources,
		false
	);

	// Now, check whether we actually got something in the kernel.
	if (!inKernel<Zp>(C.Coboundary.Matrices[2], sample, &ZZ)) RESULT = FAIL;

	resources.spinDown(&listener);
	
	return RESULT;
}
