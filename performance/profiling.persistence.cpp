
#include <ATEAMS++/ATEAMS++.h>
#include "helpers.h"

using namespace ATEAMS;
using namespace ATEAMS::topology;
using namespace std;


int main(int argc, char* argv[]) {
	string HOSTNAME = argv[1];
	int SCALE = stoi(argv[2]);
	int DIMENSION = stoi(argv[3]);
	int FIELD = stoi(argv[4]);
	int TRIALS = stoi(argv[5]);
	string STRATEGY = argv[6];

	if (DIMENSION >= 6 && SCALE >= 11) return 1;
	
	// Construct a cubical complex.
	Zp R(FIELD);
	complexes::Cubical<Zp> plex(vector<int>(DIMENSION, SCALE));
	plex.constructBoundaryMatrices(&R);
	plex.constructFullBoundaryMatrix(&R);

	// Read in the sample filtrations.
	vector<vector<int>> filtrations = filtrationData(plex.size(), SCALE, DIMENSION, TRIALS);

	// Create compute options.
	arithmetic::ComputeOptions<Zp> options;

	thread listener = options.spinUp();
	options.parallel->enabled = true;
	options.parallel->build(plex.Cells.size(), plex.size());

	// Create a bucket for storing times to completion.
	vector<int> TTC(TRIALS);
	vector<int> K;

	for (int t=0; t < TRIALS; t++) {
		K = filtrations[t];

		// ik it's lazy but whatever
		if (STRATEGY == "JIT") persistence::JIT<Zp>(&plex, K, &R, options);
		else if (STRATEGY == "twist") persistence::twist<Zp>(&plex, K, &R, options);
		else if (STRATEGY == "standard") persistence::standard<Zp>(&plex, K, &R, options);
		else if (STRATEGY == "stagger")persistence::stagger<Zp>(&plex, K, &R, options);
		else persistence::parallel<Zp>(&plex, K, &R, options);

	}

	// Spin down the listener.
	options.spinDown(&listener);

	return 0;
}