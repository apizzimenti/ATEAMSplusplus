
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
	
	// Field; thread options.
	Zp R(3);

	// RNGs and stuff.
	random_device rd;
	mt19937 RNG(rd());

	uniform_int_distribution<int> uniformValues(1, R.characteristic);
	uniform_int_distribution<int> uniformEntries(1, LENGTH);

	arithmetic::ComputeOptions options;
	thread listener = options.spinUp();

	for (int t=0; t < TRIALS; t++) {
		// Create two random vectors; check how much their indices overlap.
		int ulength = uniformEntries(RNG);
		int vlength = uniformEntries(RNG);

		SparseVector<Zp> u = randomizedVector(LENGTH, ulength, uniformValues, RNG);
		SparseVector<Zp> v = randomizedVector(LENGTH, vlength, uniformValues, RNG);
		
		if (PARALLEL) arithmetic::SparseVectorAddition(u, v, &R, options);
		else arithmetic::SparseVectorAddition(u, v, &R);
	}

	options.spinDown(&listener);

	return 0;
}
