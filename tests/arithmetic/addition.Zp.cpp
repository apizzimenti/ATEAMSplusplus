
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

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


int main(int argc, char* argv[]) {
	// Field; thread options.
	Zp R(stoi(argv[1]));
	int N = stoi(argv[2]);
	int TRIALS = stoi(argv[3]);
	int RESULT = PASS;

	// RNGs and stuff.
	random_device rd;
	mt19937 RNG(rd());

	uniform_int_distribution<int> uniformValues(1, R.characteristic-1);
	uniform_int_distribution<int> uniformEntries(1, N);

	arithmetic::ComputeOptions<Zp> options;
	thread listener = options.spinUp();


	for (int t=0; t < TRIALS; t++) {
		// Create two random vectors; check how much their indices overlap.
		int ulength = uniformEntries(RNG);
		int vlength = uniformEntries(RNG);

		SparseVector<Zp> lu = randomizedVector(N, ulength, uniformValues, RNG);
		SparseVector<Zp> ru;

		for (int i=0; i < lu.size(); i++) ru.push_back(lu(i), lu[i]);
		

		SparseVector<Zp> lv = randomizedVector(N, vlength, uniformValues, RNG);
		SparseVector<Zp> rv;

		for (int i=0; i < lv.size(); i++) rv.push_back(lv(i), lv[i]);

		SparseVector<Zp> lr = arithmetic::serialSparseVectorAddition<Zp>(lu, lv, &R);
		SparseVector<Zp> rr = arithmetic::parallelSparseVectorAddition<Zp>(ru, rv, &R, options);

		if (lr != rr) {
			RESULT = FAIL;
			print_vec_info<INDEX,typename Zp::dtype>(lu);
			print_vec_info<INDEX,typename Zp::dtype>(lv);
			print_vec_info<INDEX,typename Zp::dtype>(lr);

			cout << endl;

			print_vec_info<INDEX,typename Zp::dtype>(ru);
			print_vec_info<INDEX,typename Zp::dtype>(rv);
			print_vec_info<INDEX,typename Zp::dtype>(rr);

			cout << endl << endl << endl;
			break;
		}
	}

	options.spinDown(&listener);

	return RESULT;
}
