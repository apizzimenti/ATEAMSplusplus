
#ifndef ATEAMS_TEST_H
#define ATEAMS_TEST_H

#include <ATEAMS++/ATEAMS++.h>
#include <vector>
#include <map>

const std::vector<int> DIMENSIONS{2,3,4};

std::map<int,int> HOMOLOGICALRANK{
	{2, 2},
	{3, 3},
	{4, 6}
};

const int PASS = 0;
const int FAIL = 1;
const int ITERATIONS = 10;

inline bool inKernel(ATEAMS::ZpMatrix K, ATEAMS::ZpVector v, ATEAMS::Zp F) {
	ATEAMS::ZpVector w = SparseRREF::sparse_mat_dot_sparse_vec<ATEAMS::data_t,ATEAMS::index_t>(
		K, v, F
	);

	return w.size() < 1;
}

#endif
