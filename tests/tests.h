
#ifndef ATEAMS_TEST_H
#define ATEAMS_TEST_H

#include <ATEAMS++/ATEAMS++.h>
#include <SparseRREF/sparse_vec.h>
#include <vector>
#include <map>

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

// Defines pass/fail, default number of iterations.
const int PASS = 0;
const int FAIL = 1;
const int ITERATIONS = 10;

// Checks whether a vector is in the kernel of a matrix.
template <typename RingLike>
inline bool inKernel(ATEAMS::SparseMatrix<RingLike> K, ATEAMS::SparseVector<RingLike> v, ATEAMS::Ring* R, bool DEBUG=true) {
	ATEAMS::SparseVector<RingLike> w = SparseRREF::sparse_mat_dot_sparse_vec<typename RingLike::dtype,ATEAMS::INDEX>(
		K, v, R->ring
	);

	if (DEBUG) {
		SparseRREF::print_vec_info<ATEAMS::INDEX,typename RingLike::dtype>(v);
		SparseRREF::print_vec_info<ATEAMS::INDEX,typename RingLike::dtype>(w);
	}

	return w.size() < 1;
}

#endif
