
#ifndef ATEAMS_TEST_H
#define ATEAMS_TEST_H

#include <ATEAMS++/ATEAMS++.h>

const int PASS = 0;
const int FAIL = 1;

inline bool inKernel(ATEAMS::ZpMatrix K, ATEAMS::ZpVector v, ATEAMS::Zp F) {
	ATEAMS::ZpVector w = SparseRREF::sparse_mat_dot_sparse_vec<ATEAMS::data_t,ATEAMS::index_t>(
		K, v, F
	);

	return w.size() < 1;
}

#endif
