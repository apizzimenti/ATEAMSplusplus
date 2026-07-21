
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;


double vectorIndexOverlap(SparseVector<Zp> u, SparseVector<Zp> v) {
	set<int> uentries, ventries, intersection;

	for (int t=0; t < u.size(); t++) {
		uentries.insert((int)u(t));
		ventries.insert((int)v(t));
	}

	set_intersection(
		uentries.begin(), uentries.end(),
		ventries.begin(), ventries.end(),
		inserter(intersection, intersection.begin())
	);

	return (double)intersection.size()/(double)u.size();
}


SparseVector<Zp> randomizedVector(
	int L,
	int entries,
	uniform_int_distribution<int>& uniformValues,
	mt19937& RNG
) {
	vector<int> _random(L, 0);
	SparseVector<Zp> random;

	// Create a vector of length N*density with random entries, insert it into
	// `random`, then shuffle.
	for (int t=0; t < entries; t++) _random[t] = uniformValues(RNG);
	shuffle(_random.begin(), _random.end(), RNG);

	// Populate the SparseVector.
	for (int t=0; t < L; t++) {
		if (_random[t] > 0) {
			random.push_back((INDEX)t, (typename Zp::dtype)_random[t]);
		}
	}

	return random;
}


void writeData(string filename, string data) {
	ofstream file(filename, ios_base::app | ios_base::out);
	file << data;
}


int main(int argc, char* argv[]) {
	// Get arguments for the number of trials, the vector length.
	string HOSTNAME = argv[1];
	int TRIALS = stoi(argv[2]);
	int LENGTH = stoi(argv[3]);
	bool PARALLEL = (bool)stoi(argv[4]);

	// Column separation character.
	string sep = ",";
	
	// Field; thread options.
	Zp R(3);

	// RNGs and stuff.
	random_device rd;
	mt19937 RNG(rd());

	uniform_int_distribution<int> uniformValues(1, R.characteristic-1);
	uniform_int_distribution<int> uniformEntries(1, LENGTH);

	vector<int> lN(TRIALS);
	vector<int> rN(TRIALS);
	vector<double> OVERLAP(TRIALS);
	vector<int> TTC(TRIALS);

	arithmetic::ComputeOptions<Zp> options;
	thread listener = options.spinUp();

	for (int t=0; t < TRIALS; t++) {
		// Create two random vectors; check how much their indices overlap.
		int ulength = uniformEntries(RNG);
		int vlength = uniformEntries(RNG);

		SparseVector<Zp> u = randomizedVector(LENGTH, ulength, uniformValues, RNG);
		SparseVector<Zp> v = randomizedVector(LENGTH, vlength, uniformValues, RNG);

		long double overlap = vectorIndexOverlap(u, v);

		auto start = chrono::high_resolution_clock::now();

		if (PARALLEL) arithmetic::parallelSparseVectorAddition<Zp>(u, v, &R, options);
		else arithmetic::serialSparseVectorAddition<Zp>(u, v, &R);

		auto end = chrono::high_resolution_clock::now();

		auto duration = chrono::duration_cast<chrono::microseconds>(end-start);

		lN[t] = ulength;
		rN[t] = vlength;
		OVERLAP[t] = overlap;
		TTC[t] = duration.count();
	}

	options.spinDown(&listener);

	string csv = "";

	for (int t=0; t < TRIALS; t++) {
		csv = csv + format("{},{},{},{:.12f},{}\n", LENGTH, lN[t], rN[t], OVERLAP[t], TTC[t]);
	}

	string KIND = PARALLEL ? "parallel" : "serial";

	// APPEND to file.
	ofstream file;
	file.open(format("./performance/timing/{}.addition.{}.{}.csv", HOSTNAME, KIND, TRIALS), fstream::app);
	file << csv;
	file.close();

	return 0;
}
