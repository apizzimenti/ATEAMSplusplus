
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
	int N,
	double density,
	Ring* R
) {
	vector<int> _random(N, 0);
	SparseVector<Zp> random;

	random_device rd;
	mt19937 RNG(rd());
	uniform_int_distribution<int> intuniform(0, R->characteristic);

	// Create a vector of length N*density with random entries, insert it into
	// `random`, then shuffle.
	int entries = N*density;
	for (int t=0; t < entries; t++) _random[t] = intuniform(RNG);
	shuffle(_random.begin(), _random.end(), RNG);

	// Populate the SparseVector.
	for (int t=0; t < N; t++) {
		if (_random[t] > 0) {
			random.push_back((INDEX)t, (typename Zp::dtype)_random[t]);
		}
	}

	return random;
}


int main(int argc, char* argv[]) {
	// Get arguments for the number of trials, the vector length.
	int TRIALS = stoi(argv[1]);
	int LENGTH = stoi(argv[2]);

	// Densities are BASE^-smallestPower, in decreasing order.
	double base = sqrt(2);
	int smallestPower = 36;

	// Column separation character.
	string sep = "\t";

	// Field; thread options.
	Zp R(3);

	arithmetic::ThreadOptions options;
	options.parallelSparseAddition = false;

	thread listener = options.spinUp();

	// vector<string> columns = {"width", "density", "overlap", "mus"};

	// for (string col : columns) {
	// 	cout << format("{:<16}", col) << "\t";
	// }

	// cout << endl;

	for (int power=smallestPower; power >= 0; power--) {
		if (LENGTH*pow(base,-power) <= 5) continue;

		for (int t=0; t < TRIALS; t++) {
			// Create two random vectors; check how much their indices overlap.
			SparseVector<Zp> u = randomizedVector(LENGTH, pow(base,-power), &R);
			SparseVector<Zp> v = randomizedVector(LENGTH, pow(base,-power), &R);

			long double overlap = vectorIndexOverlap(u, v);

			auto start = chrono::high_resolution_clock::now();
			arithmetic::SparseVectorAddition(u, v, &R);
			auto end = chrono::high_resolution_clock::now();

			auto duration = chrono::duration_cast<chrono::microseconds>(end-start);
			cout << format("{:<12}", LENGTH) << sep;
			cout << format("{:<12.12f}", base) << sep;
			cout << format("{:<12}", power) << sep;
			cout << format("{:<12.12f}", pow(base,-power)) << sep;
			cout << format("{:<12.12f}", overlap) << sep;
			cout << format("{:<12}",duration.count()) << endl;
		}
	}

	options.spinDown(&listener);
	return 0;
}