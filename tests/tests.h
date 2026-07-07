
#ifndef ATEAMS_TEST_H
#define ATEAMS_TEST_H

#include <ATEAMS++/ATEAMS++.h>
#include <vector>
#include <map>

// Typical data types.
typedef ATEAMS::ff data_t;
typedef ATEAMS::SparseMatrix<data_t> ZpMatrix;
typedef ATEAMS::SparseVector<data_t> ZpVector;

const std::vector<int> DIMENSIONS{2,3,4};

std::map<int,int> HOMOLOGICALRANK{
	{2, 2},
	{3, 3},
	{4, 6}
};

const int PASS = 0;
const int FAIL = 1;
const int ITERATIONS = 10;

template <typename T>
inline bool inKernel(ATEAMS::SparseMatrix<T> K, ATEAMS::SparseVector<T> v, ATEAMS::Field F) {
	ATEAMS::SparseVector<T> w = SparseRREF::sparse_mat_dot_sparse_vec<T,ATEAMS::index_t>(
		K, v, F
	);

	return w.size() < 1;
}

#endif
