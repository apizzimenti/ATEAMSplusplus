
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

#include <SparseRREF/sparse_vec.h>

using namespace ATEAMS;
using namespace std;


bool checkReindexing(complexes::Complex<Zp>* COMPLEX, int dimension) {
	SparseMatrix<Zp> COBOUNDARY = COMPLEX->Coboundary.Full;

	// Swap two elements and verify they are reindexed correctly.
	vector<int> FILTRATION(COMPLEX->size(), 0);
	iota(FILTRATION.begin(), FILTRATION.end(), 0);

	int OFFSET = 10;
	int FIRSTINDEX = COMPLEX->Breaks[dimension][0], SECONDINDEX = COMPLEX->Breaks[dimension][0]+OFFSET;

	SparseVector<Zp> FIRST = COBOUNDARY[FIRSTINDEX];
	SparseVector<Zp> SECOND = COBOUNDARY[SECONDINDEX];

	// Find the first 3-dimensional cells with these indices in them.
	int FIRSTFACE, SECONDFACE;
	bool foundfirst = false, foundsecond = false;

	for (int t=COMPLEX->Breaks[dimension+1][0]; t < COMPLEX->Breaks[dimension+1][1]; t++) {
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
	SparseMatrix<Zp> REINDEXED = topology::reindexSparseBoundaryMatrix<Zp>(COMPLEX, FILTRATION, dimension);

	bool FIRSTREINDEXED = false, SECONDREINDEXED = false;

	for (int i=0; i < COBOUNDARY[FIRSTFACE].size(); i++) {
		if (REINDEXED[FIRSTFACE](i) == SECONDINDEX) FIRSTREINDEXED = true;
		if (REINDEXED[SECONDFACE](i) == FIRSTINDEX) SECONDREINDEXED = true;
	}

	bool SWAPPED = COBOUNDARY[FIRSTINDEX] == REINDEXED[SECONDINDEX] && COBOUNDARY[SECONDINDEX] == REINDEXED[FIRSTINDEX];

	return SWAPPED && FIRSTREINDEXED && SECONDREINDEXED;
}

template <typename RingLike>
bool checkTwistPersistence(
	complexes::Complex<RingLike>* COMPLEX,
	int dimension,
	int expectedRank,
	Ring* R,
	arithmetic::ComputeOptions<RingLike>& options
) {
	// Swap two elements and verify they are reindexed correctly.
	vector<int> FILTRATION(COMPLEX->size(), 0);
	iota(FILTRATION.begin(), FILTRATION.end(), 0);

	// We should have 1 + 4 + 6 + 4 + 1 = 16 giant components, but only 6 of the
	// desired dimension (2).
	vector<int> times = topology::twistPersistence<RingLike>(COMPLEX, FILTRATION, R, dimension, options);
	printvector<int>(times);
	return times.size() == expectedRank;
}

template <typename RingLike>
bool checkPHATPersistence(
	complexes::Complex<RingLike>* COMPLEX,
	int dimension,
	int expectedRank,
	Ring* R,
	arithmetic::ComputeOptions<RingLike>& options
) {
	// Swap two elements and verify they are reindexed correctly.
	vector<int> FILTRATION(COMPLEX->size(), 0);
	iota(FILTRATION.begin(), FILTRATION.end(), 0);

	// We should have 1 + 4 + 6 + 4 + 1 = 16 giant components, but only 6 of the
	// desired dimension (2).
	vector<int> times = topology::PHATPersistence<RingLike>(COMPLEX, FILTRATION, dimension);
	printvector<int>(times);
	return times.size() == expectedRank;
}

template <typename RingLike>
bool checkStandardPersistence(
	complexes::Complex<RingLike>* COMPLEX,
	int dimension,
	int expectedRank,
	Ring* R,
	arithmetic::ComputeOptions<RingLike>& options
) {
	// Swap two elements and verify they are reindexed correctly.
	vector<int> FILTRATION(COMPLEX->size(), 0);
	iota(FILTRATION.begin(), FILTRATION.end(), 0);

	// We should have 1 + 4 + 6 + 4 + 1 = 16 giant components, but only 6 of the
	// desired dimension (2).
	vector<int> times = topology::standardPersistence<RingLike>(COMPLEX, FILTRATION, R, dimension, options);
	printvector<int>(times);
	return times.size() == expectedRank;
}


int main(int argc, char *argv[]) {
	int RESULT = PASS;

	int FIELD = stoi(argv[1]);
	Zp ZZ(FIELD);

	// Construct arithmetic options.
	arithmetic::ComputeOptions<Zp> options;
	std::thread listener = options.spinUp();

	for (auto [dimension, rank] : HOMOLOGICALRANK) {
		// Construct a Complex.
		complexes::Cubical<Zp> COMPLEX(vector<int>(dimension, 3));

		// Construct boundary matrices.
		COMPLEX.constructBoundaryMatrices(&ZZ);
		COMPLEX.constructFlatBoundaryMatrix();
		COMPLEX.constructFullBoundaryMatrix(&ZZ);

		options.initializeParallelism(COMPLEX.Cells.size()+1);

		// Check (for fields > 2) that we're reindexing (co)boundary matrices
		// correctly. The reindexing check relies on the complex being cubical.
		if (FIELD > 2) {
			if (!checkReindexing(&COMPLEX, dimension/2)) {
				RESULT = FAIL;
				goto exit;
			}
		}

		if (!checkTwistPersistence(&COMPLEX, dimension/2, rank, &ZZ, options)) {
			RESULT = FAIL;
			goto exit;
		}

		if (!checkStandardPersistence(&COMPLEX, dimension/2, rank, &ZZ, options)) {
			RESULT = FAIL;
			goto exit;
		}
	}

	exit:
		// Now, check that we're getting the persistence correct.
		options.spinDown(&listener);
		return RESULT;
}