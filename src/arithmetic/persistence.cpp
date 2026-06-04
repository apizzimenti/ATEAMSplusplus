
#include "ATEAMS++/arithmetic/persistence.h"

#include "libraries/PHAT/compute_persistence_pairs.h"
#include "libraries/PHAT/boundary_matrix.h"
#include "libraries/PHAT/representations/default_representations.h"
#include "libraries/PHAT/algorithms/twist_reduction.h"
#include "libraries/PHAT/algorithms/standard_reduction.h"
#include "libraries/PHAT/algorithms/row_reduction.h"
#include "libraries/PHAT/algorithms/chunk_reduction.h"
#include "libraries/PHAT/algorithms/spectral_sequence_reduction.h"
#include "libraries/PHAT/algorithms/swap_twist_reduction.h"
#include "libraries/PHAT/algorithms/exhaustive_compress_reduction.h"
#include "libraries/PHAT/algorithms/lazy_retrospective_reduction.h"
#include "libraries/PHAT/helpers/dualize.h"

using namespace ATEAMS;
using namespace std;

set<int> arithmetic::PHATPersistence(arithmetic::PHATBoundaryMatrix B, vector<int> filtration) {
	set<int> r;
	return r;
}

