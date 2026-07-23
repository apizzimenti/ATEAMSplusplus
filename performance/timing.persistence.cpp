
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
	thread listener = options.spinUp(plex.Cells.size());
	options.parallel->enabled = (bool)PARALLEL;
	

	// Create a bucket for storing times to completion.
	vector<int> TTC(TRIALS);

	for (int t=0; t < TRIALS; t++) {
		// Create the filtration.
		vector<int> K = filtrate(&plex, filtration, include, DIMENSION/2);
		vector<int> times;

		auto start = chrono::high_resolution_clock::now();

		if (STRATEGY == "standard") times = topology::standardPersistence<Zp>(&plex, K, &R, DIMENSION/2, options);
		else times = topology::twistPersistence<Zp>(&plex, K, &R, DIMENSION/2, options);

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
			TTC[t],
			(int)(PARALLEL ? options.opt->pool.get_thread_count() : 1)
		);
	}

	// Spin down the listener.
	options.spinDown(&listener);

	string KIND = PARALLEL ? "parallel" : "serial";

	// APPEND to file.
	ofstream file;
	file.open(format("./performance/timing/{}.persistence.{}.{}.csv", HOSTNAME, STRATEGY, TRIALS), fstream::app);
	file << csv;
	file.close();

	return 0;
}