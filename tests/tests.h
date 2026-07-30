
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

std::map<int,int> TOTALRANKBYDIMENSION {
	{2, 4},
	{3, 8},
	{4, 16}
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


















// ################################################################################
// ## ARITHMETIC TESTS ############################################################
// ################################################################################

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
using RestrictedPersistenceAlgorithm = std::function<
		std::vector<int>
		(
			ATEAMS::complexes::Complex<RingLike>*,
			std::vector<int>&,
			ATEAMS::Ring*,
			int,
			ATEAMS::arithmetic::ComputeOptions<RingLike>&
		)
	>;

template <typename RingLike>
using FullPersistenceAlgorithm = std::function<
		std::vector<int>
		(
			ATEAMS::complexes::Complex<RingLike>*,
			std::vector<int>&,
			ATEAMS::Ring*,
			ATEAMS::arithmetic::ComputeOptions<RingLike>&
		)
	>;



template <typename RingLike>
inline bool checkRestrictedPersistence(
	ATEAMS::complexes::Complex<RingLike>* complex,
	int dimension,
	int expectedrank,
	ATEAMS::Ring* R,
	ATEAMS::arithmetic::ComputeOptions<RingLike>& options,
	std::mt19937 RNG,
	RestrictedPersistenceAlgorithm<RingLike> persistenceAlgorithm
) {
	// Initialize a filtration.
	std::vector<int> filtration(complex->size(), 0);
	std::iota(filtration.begin(), filtration.end(), 0);

	// Get the subset of things we want to shuffle.
	std::vector<int> subset(complex->Cells[dimension]);
	std::iota(subset.begin(), subset.end(), 0);

	// Shuffle them, and insert into the filtration.
	std::shuffle(subset.begin(), subset.end(), RNG);
	int offset = (dimension > 0) ? complex->Offsets[dimension-1] : 0;

	for (int t=complex->Offsets[dimension-1]; t < complex->Offsets[dimension]; t++) {
		filtration[t] = subset[t-offset]+offset;
	}

	// Check whether the rank is correct.
	std::vector<int> times = persistenceAlgorithm(complex, filtration, R, dimension, options);
	printvector<int>(times);
	return times.size() == expectedrank;
}


template <typename RingLike>
inline bool checkFullPersistence(
	ATEAMS::complexes::Complex<RingLike>* complex,
	int dimension,
	int expectedTotalRank,
	ATEAMS::Ring* R,
	ATEAMS::arithmetic::ComputeOptions<RingLike>& options,
	std::mt19937 RNG,
	FullPersistenceAlgorithm<RingLike> persistenceAlgorithm
) {
	// Initialize a filtration.
	std::vector<int> filtration(complex->size(), 0);
	std::iota(filtration.begin(), filtration.end(), 0);

	// Get the subset of things we want to shuffle.
	std::vector<int> subset(complex->Cells[dimension]);
	std::iota(subset.begin(), subset.end(), 0);

	// Shuffle them, and insert into the filtration.
	std::shuffle(subset.begin(), subset.end(), RNG);
	int offset = (dimension > 0) ? complex->Offsets[dimension-1] : 0;

	for (int t=complex->Offsets[dimension-1]; t < complex->Offsets[dimension]; t++) {
		filtration[t] = subset[t-offset]+offset;
	}

	// Check whether the rank is correct.
	std::vector<int> times = persistenceAlgorithm(complex, filtration, R, options);
	printvector<int>(times);
	return times.size() == expectedTotalRank;
}


template <typename RingLike>
inline bool checkSingleReindexing(
	ATEAMS::complexes::Complex<RingLike>* complex,
	int dimension,
	ATEAMS::arithmetic::ComputeOptions<RingLike>& options
) {
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

	for (int t=0; t < complex->size(); t++) {
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
	ATEAMS::SparseMatrix<RingLike> FullReindexed = ATEAMS::topology::persistence::policies::singleReindexingPolicy<RingLike>(complex, filtration, options, dimension);

	bool firstReindexed = false, secondReindexed = false;

	for (int i=0; i < Full[firstFace].size(); i++) {
		if (FullReindexed[firstFace](i) == secondIndex) firstReindexed = true;
		if (FullReindexed[secondFace](i) == firstIndex) secondReindexed = true;
	}

	bool SWAPPED = Full[firstIndex] == FullReindexed[secondIndex] && Full[secondIndex] == FullReindexed[firstIndex];

	return SWAPPED && firstReindexed && secondReindexed;
}


template <typename RingLike>
inline bool checkFullReindexing(
	ATEAMS::complexes::Complex<RingLike>* complex,
	ATEAMS::arithmetic::ComputeOptions<RingLike>& options
) {
	ATEAMS::SparseMatrix<RingLike> Full = complex->Coboundary.Full;

	// Swap two elements and verify they are reindexed correctly.
	std::vector<int> filtration(complex->size(), 0);
	std::iota(filtration.begin(), filtration.end(), 0);

	int dimension = 1;

	int offset = 10;
	int firstIndex = complex->Breaks[dimension][0], secondIndex = complex->Breaks[dimension][0]+offset;

	ATEAMS::SparseVector<RingLike> first = Full[firstIndex];
	ATEAMS::SparseVector<RingLike> second = Full[secondIndex];

	// Find the first 3-dimensional cells with these indices in them.
	int firstFace, secondFace;
	bool foundfirst = false, foundsecond = false;

	for (int t=0; t < complex->size(); t++) {
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
	ATEAMS::SparseMatrix<RingLike> FullReindexed = ATEAMS::topology::persistence::policies::singleReindexingPolicy<RingLike>(complex, filtration, options, dimension);

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
	RestrictedPersistenceAlgorithm<RingLike> restrictedPersistenceAlgorithm,
	FullPersistenceAlgorithm<RingLike> fullPersistenceAlgorithm,
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
		
		// Make sure we've set the parallel computing options correctly.
		options.parallel->enabled = true;
		options.parallel->build(complex.size(), complex.Cells.size());
		options.serial->build(complex.size());

		// Check whether we're re-indexing properly.
		if (!checkSingleReindexing<RingLike>(&complex, dimension/2, options)) {
			RESULT = FAIL;
			goto EXIT;
		}

		// Check whether we're persisting properly.
		for (int t=0; t < 256; t++) {
			if (!checkRestrictedPersistence<RingLike>(&complex, dimension/2, rank, &R, options, RNG, restrictedPersistenceAlgorithm)) {
				RESULT = FAIL;
				goto EXIT;
			}

			if (!checkFullPersistence<RingLike>(&complex, dimension/2, TOTALRANKBYDIMENSION[dimension], &R, options, RNG, fullPersistenceAlgorithm)) {
				RESULT = FAIL;
				goto EXIT;
			}
		}
	}

	EXIT:
		options.spinDown(&listener);
		return RESULT;
}

#endif
