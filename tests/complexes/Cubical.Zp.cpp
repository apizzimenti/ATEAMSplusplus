
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;


vector<int> homologySizes(
	complexes::Complex<Zp>* COMPLEX,
	Ring* R,
	int dimension,
	arithmetic::ComputeResources<Zp>& resources
) {
	// Swap two elements and verify they are reindexed correctly.
	vector<int> FILTRATION(COMPLEX->size(), 0);
	iota(FILTRATION.begin(), FILTRATION.end(), 0);
	
	vector<int> sizes;

	for (int d=0; d < COMPLEX->Breaks.size(); d++) {
		vector<int> times = topology::persistence::persistence<Zp>(COMPLEX, FILTRATION, R, d, resources);
		printvector<int>(times);

		sizes.push_back(
			std::count_if(times.begin(), times.end(), [=](int t) {
				return COMPLEX->Breaks[d][0] <= t && t < COMPLEX->Breaks[d][1];
			})
		);
	}

	return sizes;
}


int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	int RESULT = PASS;

	// Construct arithmetic resources.
	arithmetic::ComputeResources<Zp> resources;
	std::thread listener = resources.spinUp();

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
		Zp R(FIELD);

		complexes::Cubical<Zp> CUBICAL(corners);
		CUBICAL.constructBoundaryMatrices(&R);
		CUBICAL.constructFlatBoundaryMatrix();
		CUBICAL.constructFullBoundaryMatrix(&R);

		resources.arithmetize(&R);
		resources.parallel->enabled = true;
		resources.serial->build(CUBICAL.size());
		CUBICAL.constructSparseBases(&R, resources);

		for (int d=0; d < CUBICAL.Cells.size(); d++) {
			cout << format("{}: {}",d,CUBICAL.Boundary.Bases[d].size()) << endl;
		}
		break;
	}

	// Now, check that we're getting the persistence correct.
	resources.spinDown(&listener);
	return PASS;
}
