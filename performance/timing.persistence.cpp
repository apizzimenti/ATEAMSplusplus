
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

	// Construct a cubical complex and the ingredients for a filtration.
	Zp R(FIELD);
	complexes::Cubical<Zp> plex(vector<int>(DIMENSION, SCALE));
	plex.constructBoundaryMatrices(&R);
	plex.constructFullBoundaryMatrix(&R);

	// Read in the sample filtrations.
	vector<vector<int>> filtrations = filtrationData(plex.size(), SCALE, DIMENSION, TRIALS);
	
	// Create compute resources.
	arithmetic::ComputeResources<Zp> resources;

	thread listener = resources.spinUp();
	resources.arithmetize(&R);
	resources.parallel->enabled = true;
	resources.parallel->build(plex.size(), plex.Cells.size());
	resources.serial->build(plex.size());


	// Create a bucket for storing times to completion.
	vector<int> TTC(TRIALS);

	for (int t=0; t < TRIALS; t++) {
		// Create the filtration.
		vector<int> K = filtrations[t];

		auto start = chrono::high_resolution_clock::now();

		// ik it's lazy but whatever
		if (STRATEGY == "JIT") persistence::JIT<Zp>(&plex, K, &R, resources);
		else if (STRATEGY == "twist") persistence::twist<Zp>(&plex, K, &R, resources);
		else if (STRATEGY == "standard") persistence::standard<Zp>(&plex, K, &R, resources);
		else if (STRATEGY == "stagger") persistence::stagger<Zp>(&plex, K, &R, resources);
		else if (STRATEGY == "split") persistence::split<Zp>(&plex, K, &R, resources);
		else persistence::parallel<Zp>(&plex, K, &R, resources);

		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::microseconds>(end-start);

		// Write data to storage.
		TTC[t] = duration.count();
	}

	// Write data to file.
	string csv = "";

	for (int t=0; t < TRIALS; t++) {
		csv = csv + format(
			"{},{},{},{}\n",
			SCALE,
			DIMENSION,
			FIELD,
			TTC[t]
		);
	}

	// Spin down the listener.
	resources.spinDown(&listener);

	// APPEND to file.
	ofstream file;
	file.open(format("./performance/timing/{}.persistence.{}.{}.csv", HOSTNAME, STRATEGY, TRIALS), fstream::app);
	file << csv;
	file.close();

	return 0;
}