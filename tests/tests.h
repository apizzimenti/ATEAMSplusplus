
#ifndef ATEAMS_TEST_H
#define ATEAMS_TEST_H

#include <ATEAMS++/ATEAMS++.h>
#include <vector>
#include <map>

// Typical data types.
typedef ATEAMS::ff FINITE;
typedef ATEAMS::rational RATIONAL;

// Dimensions over which we're testing.
const std::vector<int> DIMENSIONS{2,3,4};

// Correct ranks of homology groups.
std::map<int,int> HOMOLOGICALRANK{
	{2, 2},
	{3, 3},
	{4, 6}
};

// Defines pass/fail, default number of iterations.
const int PASS = 0;
const int FAIL = 1;
const int ITERATIONS = 10;

// Checks whether a vector is in the kernel of a matrix.
template <typename T>
inline bool inKernel(ATEAMS::SparseMatrix<T> K, ATEAMS::SparseVector<T> v, ATEAMS::Field F) {
	ATEAMS::SparseVector<T> w = SparseRREF::sparse_mat_dot_sparse_vec<T,ATEAMS::index_t>(
		K, v, F
	);

	return w.size() < 1;
}

#endif
