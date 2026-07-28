
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace ATEAMS::topology;
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
	int TRIALS = stoi(argv[3]);

	// Construct a cubical complex and the ingredients for a filtration.
	Zp R(2);
	complexes::Cubical<Zp> plex(vector<int>(DIMENSION, SCALE));
	plex.constructBoundaryMatrices(&R);
	plex.constructFullBoundaryMatrix(&R);

	// Create reusable indices for filtrations.
	vector<int> filtration(plex.size(), 0);
	iota(begin(filtration), end(filtration), 0);

	vector<vector<int>> includes(plex.Cells.size());

	for (int d=0; d < plex.Cells.size(); d++) {
		vector<int> include(plex.Cells[d], 0);
		iota(begin(include), end(include), 0);

		includes[d] = include;
	}

	// Generate a bunch of filtrations.
	vector<vector<int>> filtrations(TRIALS);

	#pragma omp parallel for shared(filtrations) firstprivate(includes)
	for (int t=0; t < TRIALS; t++) {
		// Create reusable indices for filtrations.
		vector<int> filtration(plex.size(), 0);
		iota(begin(filtration), end(filtration), 0);

		for (int d=0; d < plex.Breaks.size(); d++) {
			filtration = filtrate(&plex, filtration, includes[d], d);
		}

		filtrations[t] = filtration;
	}

	// Write data to file.
	string csv = "";

	for (int t=0; t < TRIALS; t++) {
		for (int i=0; i < filtrations[t].size()-1; i++) {
			csv += format("{},", filtrations[t][i]);
		}
		csv += format("{}\n",filtrations[t][filtrations[t].size()-1]);
	}

	// APPEND to file.
	ofstream file;
	file.open(format("./performance/samples/{}.{}.csv", SCALE, DIMENSION));
	file << csv;
	file.close();

	return 0;
}