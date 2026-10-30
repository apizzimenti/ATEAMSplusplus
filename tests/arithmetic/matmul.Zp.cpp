
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

int main(int argc, char* argv[]) {
	// Field; thread resources.
	Zp R(stoi(argv[1]));
	int M = stoi(argv[2]);
	int N = stoi(argv[3]);
	int TRIALS = stoi(argv[4]);
	int RESULT = PASS;

	// RNGs and stuff.
	random_device rd;
	mt19937 RNG(rd());

	uniform_int_distribution<int> uniformValues(0, R.characteristic-1);

	// Construct arithmetic resources.
	ATEAMS::arithmetic::ComputeResources<Zp> resources;
	resources.arithmetize(&R);
	resources.parallel->enabled = true;
	
	std::thread listener = resources.spinUp();

	for (int t=0; t < TRIALS; t++) {
		// Create a random vector and a random matrix.
		SparseVector<Zp> x = randomizedVector(N, uniformValues, RNG);
		SparseMatrix<Zp> A = randomizedMatrix(M, N, uniformValues, RNG);

		SparseVector<Zp> sy = arithmetic::SparseRightMultiplication<Zp>(A, x, &R);
		SparseVector<Zp> py = arithmetic::OpenMPSparseRightMultiplication<Zp>(A, x, &R, resources);
		SparseVector<Zp> by = arithmetic::SparseRightMultiplication<Zp>(A, x, &R, resources);

		if (sy != py || sy != by) {
			RESULT = FAIL;
			print_vec_info<INDEX,typename Zp::dtype>(sy);
			cerr << endl;
			print_vec_info<INDEX,typename Zp::dtype>(py);
			cerr << endl;
			print_vec_info<INDEX,typename Zp::dtype>(by);
			cerr << endl;
			break;
		}
	}

	resources.spinDown(&listener);

	return RESULT;
}
