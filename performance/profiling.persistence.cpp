
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

vector<int> filtrate(
	complexes::Complex<Zp>* complex,
	vector<int>& filtration,
	vector<int>& include,
	int DIMENSION
) {
	random_device rd;
	mt19937 RNG(rd());

	shuffle(include.begin(), include.end(), RNG);

	int start = complex->Breaks[DIMENSION][0];
	int stop = complex->Breaks[DIMENSION][1];
	int offset = (DIMENSION > 0) ? complex->Offsets[DIMENSION-1] : 0;

	for (int j=0; j < include.size(); j++) filtration[j+start] = include[j]+start;

	return filtration;
}


int main(int argc, char* argv[]) {
	string HOSTNAME = argv[1];
	int SCALE = stoi(argv[2]);
	int DIMENSION = stoi(argv[3]);
	int FIELD = stoi(argv[4]);
	int TRIALS = stoi(argv[5]);
	string STRATEGY = argv[6];
	int PARALLEL = stoi(argv[7]);

	// Construct a cubical complex and the ingredients for a filtration.
	Zp R(FIELD);
	complexes::Cubical<Zp> plex(vector<int>(DIMENSION, SCALE));
	plex.constructBoundaryMatrices(&R);
	plex.constructFullBoundaryMatrix(&R);

	// Create reusable indices for filtrations.
	vector<int> filtration(plex.size(), 0);
	iota(begin(filtration), end(filtration), 0);

	vector<int> include(plex.Cells[DIMENSION/2], 0);
	iota(begin(include), end(include), 0);

	// Create compute options.
	arithmetic::ComputeOptions<Zp> options;
	thread listener = options.spinUp();

	// Make sure we've set the parallel computing options correctly.
	options.parallel->enabled = (bool)PARALLEL;
	if ((bool)PARALLEL) options.parallel->build(plex.Cells.size(), plex.size());

	vector<int> K = filtrate(&plex, filtration, include, DIMENSION/2);

	for (int t=0; t < TRIALS; t++) {
		if (STRATEGY == "standard")topology::standardPersistence<Zp>(&plex, K, &R, DIMENSION/2, options);
		else topology::twistPersistence<Zp>(&plex, K, &R, DIMENSION/2, options);
	}

	return 0;
}