
#ifndef ATEAMS_ARITHMETIC_PERSISTENCE_H
#define ATEAMS_ARITHMETIC_PERSISTENCE_H

#include <phat/compute_persistence_pairs.h>
#include <phat/boundary_matrix.h>
#include <phat/representations/default_representations.h>

#include <vector>
#include <set>

namespace ATEAMS::arithmetic {
	typedef std::vector<phat::index> PHATColumn;
	typedef phat::boundary_matrix<phat::bit_tree_pivot_column> PHATBoundaryMatrix;
	typedef phat::persistence_pairs PHATPairs;
	typedef phat::twist_reduction PHATTwist;

	std::set<int> PHATPersistence(PHATBoundaryMatrix B, std::vector<int> filtration);
}

#endif

