
#ifndef ATEAMS_TEST_H
#define ATEAMS_TEST_H

#include <ATEAMS++/ATEAMS++.h>
#include <SparseRREF/sparse_vec.h>
#include <vector>
#include <map>
#include <string>

// Typical data types.
typedef ATEAMS::FINITE FINITE;
typedef ATEAMS::RATIONAL RATIONAL;

// Dimensions over which we're testing.
const std::vector<int> DIMENSIONS{2,3,4};

// Correct ranks of homology groups.
std::map<int,int> HOMOLOGICALRANK{
	{2, 2},
	{3, 3},
	{4, 6}
};

// When to stop invading.
std::map<int,std::vector<int>> STOPINVADING{
	{2,{1}},
	{3,{1,2}},
	{4,{3,4}}
};

// Defines pass/fail, default number of iterations.
const int PASS = 0;
const int FAIL = 1;
const int ITERATIONS = 10;

// Checks whether a vector is in the kernel of a matrix.
template <typename RingLike>
inline bool inKernel(ATEAMS::SparseMatrix<RingLike> K, ATEAMS::SparseVector<RingLike> v, ATEAMS::Ring* R, bool DEBUG=true) {

	ATEAMS::SparseVector<RingLike> w = ATEAMS::arithmetic::SparseRightMultiplication<RingLike>(
		K, v, R
	);

	if (DEBUG) {
		SparseRREF::print_vec_info<ATEAMS::INDEX,typename RingLike::dtype>(v);
		SparseRREF::print_vec_info<ATEAMS::INDEX,typename RingLike::dtype>(w);
	}

	return w.size() < 1;
}


// ################################################################################
// ## PERSISTENCE TESTS ###########################################################
// ################################################################################

template <typename RingLike>
inline bool checkPersistence(
	ATEAMS::complexes::Complex<RingLike>* complex,
	int dimension,
	int expectedrank,
	ATEAMS::Ring* R,
	ATEAMS::arithmetic::ComputeOptions<RingLike>& options,
	std::mt19937 RNG,
	std::function<
		std::vector<int>
		(
			ATEAMS::complexes::Complex<RingLike>*,
			std::vector<int>&,
			ATEAMS::Ring*,
			int,
			ATEAMS::arithmetic::ComputeOptions<RingLike>&
		)
	> persistenceAlgorithm
) {
	std::vector<int> filtration(complex->size(), 0);
	std::iota(filtration.begin(), filtration.end(), 0);

	std::vector<int> subset(complex->Cells[dimension]);
	std::iota(subset.begin(), subset.end(), 0);

	std::shuffle(subset.begin(), subset.end(), RNG);
	int offset = (dimension > 0) ? complex->Offsets[dimension-1] : 0;

	for (int t=complex->Offsets[dimension-1]; t < complex->Offsets[dimension]; t++) {
		filtration[t] = subset[t-offset]+offset;
	}

	std::vector<int> times = persistenceAlgorithm(complex, filtration, R, dimension, options);
	return times.size() == expectedrank;
}

template <typename RingLike>
inline bool checkReindexing(ATEAMS::complexes::Complex<RingLike>* complex, int dimension) {
	ATEAMS::SparseMatrix<RingLike> Full = complex->Coboundary.Full;

	// Swap two elements and verify they are reindexed correctly.
	std::vector<int> filtration(complex->size(), 0);
	std::iota(filtration.begin(), filtration.end(), 0);

	int offset = 10;
	int firstIndex = complex->Breaks[dimension][0], secondIndex = complex->Breaks[dimension][0]+offset;

	ATEAMS::SparseVector<RingLike> first = Full[firstIndex];
	ATEAMS::SparseVector<RingLike> second = Full[secondIndex];

	// Find the first 3-dimensional cells with these indices in them.
	int firstFace, secondFace;
	bool foundfirst = false, foundsecond = false;

	for (int t=complex->Breaks[dimension+1][0]; t < complex->Breaks[dimension+1][1]; t++) {
		for (int i=0; i < Full[t].size(); i++) {
			if (Full[t](i) == firstIndex && !foundfirst) {
				firstFace = t;
				foundfirst = true;
			}

			if (Full[t](i) == secondIndex && !foundsecond) {
				secondFace = t;
				foundsecond = true;
			}

			if (foundfirst && foundsecond) break;
		}
	}

	// Swap, then reindex.
	filtration[firstIndex] = secondIndex;
	filtration[secondIndex] = firstIndex;
	ATEAMS::SparseMatrix<RingLike> FullReindexed = ATEAMS::topology::reindexSparseBoundaryMatrix<RingLike>(complex, filtration, dimension);

	bool firstReindexed = false, secondReindexed = false;

	for (int i=0; i < Full[firstFace].size(); i++) {
		if (FullReindexed[firstFace](i) == secondIndex) firstReindexed = true;
		if (FullReindexed[secondFace](i) == firstIndex) secondReindexed = true;
	}

	bool SWAPPED = Full[firstIndex] == FullReindexed[secondIndex] && Full[secondIndex] == FullReindexed[firstIndex];

	return SWAPPED && firstReindexed && secondReindexed;
}


template <typename RingLike>
inline int persistenceDispatcher(
	int argc,
	char *argv[],
	std::function<
		std::vector<int>
		(
			ATEAMS::complexes::Complex<RingLike>*,
			std::vector<int>&,
			ATEAMS::Ring*,
			int,
			ATEAMS::arithmetic::ComputeOptions<RingLike>&
		)
	> persistenceAlgorithm,
	bool parallel=true
) {
	int RESULT = PASS;
	int FIELD = std::stoi(argv[1]);
	RingLike R(FIELD);

	// Construct arithmetic options.
	ATEAMS::arithmetic::ComputeOptions<RingLike> options;
	std::thread listener = options.spinUp();

	// Create the RNG.
	std::random_device rd;
	std::mt19937 RNG(rd());
	
	for (auto [dimension, rank] : HOMOLOGICALRANK) {
		// Construct a Complex.
		ATEAMS::complexes::Cubical<RingLike> complex(std::vector<int>(dimension, 3));

		// Construct boundary matrices.
		complex.constructBoundaryMatrices(&R);
		complex.constructFlatBoundaryMatrix();
		complex.constructFullBoundaryMatrix(&R);

		options.initializeParallelism(complex.Cells.size()+1);
		options.parallel->enabled = parallel;

		// Check whether we're re-indexing properly.
		if (!checkReindexing<RingLike>(&complex, dimension/2)) {
			RESULT = FAIL;
			goto EXIT;
		}

		// Check whether we're persisting properly.
		for (int t=0; t < 100; t++) {
			if (!checkPersistence<RingLike>(&complex, dimension/2, rank, &R, options, RNG, persistenceAlgorithm)) {
				RESULT = FAIL;
				goto EXIT;
			}
		}
	}

	EXIT:
		// Now, check that we're getting the persistence correct.
		options.spinDown(&listener);
		return RESULT;
}

#endif
