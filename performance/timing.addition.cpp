
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;


double vectorIndexOverlap(SparseVector<Zp> u, SparseVector<Zp> v) {
	vector<int> uentries(u.size()), ventries(v.size()), intersection;

	for (int t=0; t < u.size(); t++) uentries[t] = (int)u(t);
	for (int t=0; t < v.size(); t++) ventries[t] = (int)v(t);

	set_intersection(
		uentries.begin(), uentries.end(),
		ventries.begin(), ventries.end(),
		back_inserter(intersection)
	);

	return (double)intersection.size()/(double)u.size();
}


vector<SparseVector<Zp>> randomizedVectors(
	int L,
	int entries,
	mt19937& RNG
) {
	vector<int> _random(L, 0);
	SparseVector<Zp> lrandom;
	SparseVector<Zp> rrandom;

	// Create a vector of length N*density with all-1 entries.
	for (int t=0; t < entries; t++) _random[t] = 1;

	// Populate the first SparseVector.
	shuffle(_random.begin(), _random.end(), RNG);
	for (int t=0; t < L; t++) {
		if (_random[t] > 0) {
			lrandom.push_back((INDEX)t, (typename Zp::dtype)_random[t]);
		}
	}

	// Populate the second SparseVector.
	shuffle(_random.begin(), _random.end(), RNG);
	for (int t=0; t < L; t++) {
		if (_random[t] > 0) {
			rrandom.push_back((INDEX)t, (typename Zp::dtype)_random[t]);
		}
	}

	vector<SparseVector<Zp>> cont = {lrandom, rrandom};
	return cont;
}

int main(int argc, char* argv[]) {
	// Get arguments for the number of trials, the vector length.
	string HOSTNAME = argv[1];
	int TRIALS = stoi(argv[2]);
	int LENGTH = stoi(argv[3]);
	bool PARALLEL = (bool)stoi(argv[4]);
	int STRATEGY = stoi(argv[5]);

	// Column separation character.
	string sep = ",";
	
	// Field; thread options.
	Zp R(3);

	// RNGs and stuff.
	random_device rd;
	mt19937 RNG(rd());

	uniform_int_distribution<int> uniformLength(1, LENGTH);

	vector<int> l(TRIALS);
	vector<double> OVERLAP(TRIALS);
	vector<int> TTC(TRIALS);

	arithmetic::ComputeOptions<Zp> options;
	thread listener = options.spinUp(STRATEGY);

	for (int t=0; t < TRIALS; t++) {
		// Create two random vectors; check how much their indices overlap.
		int length = uniformLength(RNG);
		vector<SparseVector<Zp>> randoms = randomizedVectors(LENGTH, length, RNG);

		long double overlap = vectorIndexOverlap(randoms[0], randoms[1]);

		auto start = chrono::high_resolution_clock::now();

		if (PARALLEL) arithmetic::parallelSparseVectorAddition<Zp>(randoms[0], randoms[1], &R, options);
		else arithmetic::serialSparseVectorAddition<Zp>(randoms[0], randoms[1], &R);

		auto end = chrono::high_resolution_clock::now();

		auto duration = chrono::duration_cast<chrono::microseconds>(end-start);

		l[t] = length;
		OVERLAP[t] = overlap;
		TTC[t] = duration.count();
	}

	options.spinDown(&listener);

	string csv = "";

	for (int t=0; t < TRIALS; t++) {
		csv = csv + format(
			"{},{},{:.12f},{},{}\n",
			LENGTH,
			l[t],
			OVERLAP[t],
			TTC[t],
			PARALLEL ? options.parallel->threads : 1
		);
	}

	string KIND = PARALLEL ? "parallel" : "serial";

	// APPEND to file.
	ofstream file;
	file.open(format("./performance/timing/{}.addition.{}.{}.csv", HOSTNAME, KIND, TRIALS), fstream::app);
	file << csv;
	file.close();

	return 0;
}
