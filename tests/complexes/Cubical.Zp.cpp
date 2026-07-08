
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;


vector<int> homologySizes(complexes::Complex<FINITE>* COMPLEX, int dimension) {
	// Swap two elements and verify they are reindexed correctly.
	vector<int> FILTRATION(COMPLEX->size(), 0);
	iota(FILTRATION.begin(), FILTRATION.end(), 0);

	// Compute 
	vector<int> times = arithmetic::PHATPersistence<FINITE>(COMPLEX, FILTRATION, dimension/2);
	vector<int> sizes;

	for (int d=0; d < COMPLEX->breaks.size(); d++) {
		sizes.push_back(
			std::count_if(times.begin(), times.end(), [=](int t) {
				return COMPLEX->breaks[d][0] <= t && t < COMPLEX->breaks[d][1];
			})
		);
	}

	return sizes;
}


int main() {
	// Construct arithmetic options.
	arithmetic::ThreadOptions options;
	std::thread listener = options.spinUp();

	// Construct Cubical complexes of varying dimensions/boundary conditions,
	// verifying whether they're constructed correctly. We perform two checks:
	// 
	//	1. whether each d-cell belongs to the appropriate number of (d+1)-cells;
	//	2. whether the ranks of the homology groups of the space are correct.
	//
	// We only need to check these for periodic boundary conditions, since their
	// correctness implies the correctness of the complex with free boundary
	// conditions.
	vector<int> dimensions{2,3,4};

	map<int,vector<int>> homologies{
		{2, {1,2,1}},
		{3, {1,3,3,1}},
		{4, {1,4,6,4,1}}
	};

	for (int dimension : dimensions) {
		vector<int> corners(dimension, 3);
		complexes::Cubical CUBICAL(corners);
		CUBICAL.constructFlatBoundaryMatrix();

		vector<int> sizes = homologySizes(&CUBICAL, dimension);

		if (sizes != homologies[dimension]) {
			options.spinDown(&listener);
			return FAIL;
		}
	}

	// Now, check that we're getting the persistence correct.
	options.spinDown(&listener);
	return PASS;
}
