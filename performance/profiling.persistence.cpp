
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
	int SCALE = stoi(argv[1]);
	int DIMENSION = stoi(argv[2]);
	int FIELD = stoi(argv[3]);
	int ATTEMPTS = stoi(argv[4]);
	int PARALLEL = stoi(argv[5]);

	// Construct a cubical complex and the ingredients for a filtration.
	Zp R(FIELD);
	complexes::Cubical<Zp> plex(vector<int>(DIMENSION*2, SCALE));
	
	plex.constructBoundaryMatrices(&R);
	plex.constructFlatBoundaryMatrix();
	plex.constructFullBoundaryMatrix(&R);

	// Create reusable indices for filtrations.
	vector<int> filtration(plex.size());
	iota(begin(filtration), end(filtration), 0);

	vector<int> include(plex.Cells[DIMENSION]);
	iota(begin(include), end(include), 0);

	arithmetic::ThreadOptions options;
	options.parallelSparseAddition = (bool)PARALLEL;
	thread listener = options.spinUp();

	for (int t=0; t < ATTEMPTS; t++) {
		// Create the filtration.
		vector<int> K = filtrate(&plex, filtration, include, DIMENSION);
		topology::persistence<Zp>(&plex, K, &R, DIMENSION, options);
	}

	options.spinDown(&listener);

	return 0;
}