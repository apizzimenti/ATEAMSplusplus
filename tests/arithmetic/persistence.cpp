
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

#include <SparseRREF/sparse_vec.h>

using namespace ATEAMS;
using namespace std;


bool checkReindexing(complexes::Complex<data_t>* COMPLEX) {
	ZpMatrix COBOUNDARY = COMPLEX->Coboundary.Full;

	// Swap two elements and verify they are reindexed correctly.
	vector<int> FILTRATION(COMPLEX->size(), 0);
	iota(FILTRATION.begin(), FILTRATION.end(), 0);

	int OFFSET = 100;
	int FIRSTINDEX = COMPLEX->breaks[2][0], SECONDINDEX = COMPLEX->breaks[2][0]+OFFSET;

	ZpVector FIRST = COBOUNDARY[FIRSTINDEX];
	ZpVector SECOND = COBOUNDARY[SECONDINDEX];

	// Find the first 3-dimensional cells with these indices in them.
	int FIRSTFACE, SECONDFACE;
	bool foundfirst = false, foundsecond = false;

	for (int t=COMPLEX->breaks[3][0]; t < COMPLEX->breaks[3][1]; t++) {
		for (int i=0; i < COBOUNDARY[t].size(); i++) {
			if (COBOUNDARY[t](i) == FIRSTINDEX && !foundfirst) {
				FIRSTFACE = t;
				foundfirst = true;
			}

			if (COBOUNDARY[t](i) == SECONDINDEX && !foundsecond) {
				SECONDFACE = t;
				foundsecond = true;
			}

			if (foundfirst && foundsecond) break;
		}
	}

	// Swap, then reindex.
	FILTRATION[FIRSTINDEX] = SECONDINDEX;
	FILTRATION[SECONDINDEX] = FIRSTINDEX;
	ZpMatrix REINDEXED = arithmetic::reindexSparseBoundaryMatrix<data_t>(COMPLEX, FILTRATION, 2);

	bool FIRSTREINDEXED = false, SECONDREINDEXED = false;

	for (int i=0; i < COBOUNDARY[FIRSTFACE].size(); i++) {
		if (REINDEXED[FIRSTFACE](i) == SECONDINDEX) FIRSTREINDEXED = true;
		if (REINDEXED[SECONDFACE](i) == FIRSTINDEX) SECONDREINDEXED = true;
	}

	bool SWAPPED = COBOUNDARY[FIRSTINDEX] == REINDEXED[SECONDINDEX] && COBOUNDARY[SECONDINDEX] == REINDEXED[FIRSTINDEX];

	return SWAPPED && FIRSTREINDEXED && SECONDREINDEXED;
}


bool checkPHATpersistence(complexes::Complex<data_t>* COMPLEX) {
	// Swap two elements and verify they are reindexed correctly.
	vector<int> FILTRATION(COMPLEX->size(), 0);
	iota(FILTRATION.begin(), FILTRATION.end(), 0);

	// We should have 1 + 4 + 6 + 4 + 1 = 16 giant components.
	vector<int> times = arithmetic::PHATPersistence<data_t>(COMPLEX, FILTRATION, 2);

	return times.size() == 16;
}

bool checkTwistPersistence(complexes::Complex<data_t>* COMPLEX, Field F) {
	// Swap two elements and verify they are reindexed correctly.
	vector<int> FILTRATION(COMPLEX->size(), 0);
	iota(FILTRATION.begin(), FILTRATION.end(), 0);

	// We should have 1 + 4 + 6 + 4 + 1 = 16 giant components, but only 6 of the
	// desired dimension (2).
	vector<int> times = arithmetic::TwistPersistence<data_t>(COMPLEX, FILTRATION, F, 2);

	return times.size() == 6;
}


int main(int argc, char *argv[]) {
	int FIELD = stoi(argv[1]);
	Field F(SparseRREF::FIELD_Fp, FIELD);

	// Construct a Complex.
	complexes::Cubical<data_t> COMPLEX({3,3,3,3});

	// Construct boundary matrices.
	COMPLEX.constructBoundaryMatrices(F);
	COMPLEX.constructFlatBoundaryMatrix();
	COMPLEX.constructFullBoundaryMatrix(F);

	// Construct arithmetic options.
	arithmetic::ThreadOptions options;
	std::thread listener = options.spinUp();

	// Check (for fields > 2) that we're reindexing (co)boundary matrices
	// correctly. The reindexing check relies on the complex being cubical.
	if (FIELD > 2) {
		if (!checkReindexing(&COMPLEX)) {
			options.spinDown(&listener);
			return FAIL;
		}
	}

	if (FIELD > 2) {
		if (!checkTwistPersistence(&COMPLEX, F)) {
			options.spinDown(&listener);
			return FAIL;
		}
	} else {
		if (!checkPHATpersistence(&COMPLEX)) {
			options.spinDown(&listener);
			return FAIL;
		}
	}

	// Now, check that we're getting the persistence correct.
	options.spinDown(&listener);
	return PASS;
}