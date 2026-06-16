
#include "ATEAMS++/common.h"
#include "ATEAMS++/util.h"
#include "ATEAMS++/complexes/Complex.h"
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

vector<int> arithmetic::PHATPersistence(ATEAMS::complexes::Complex* complex, vector<int> filtration, int dimension) {
	// The filtration specifies the order in which we add the cells of all
	// dimensions. Create a map that specifies to which position each cell was
	// moved. For example, if the filtration has
	//
	//				[... 12, 9, 10, 19, ...]
	//
	// that was originally
	//				[... 9, 10, 11, 12, ...]
	//
	// then the mapping should have entries
	//
	//				{... 12: 9, 9: 10, 10: 11, 19: 12, ...}
	map<int,int> remapping;
	for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

	PHATBoundaryMatrix boundary;
	PHATColumn column;
	int faces;

	int startDimension = complex->breaks[dimension][0];
	int stopDimension = complex->breaks[dimension][1];

	boundary.set_num_cols(complex->Boundary.Flat.size());

	for (int t=0; t < complex->Boundary.Flat.size(); t++) {
		faces = complex->Boundary.Flat[t].size();

		// Check whether we're in the dimension of the cells on which we're
		// percolating; if so, add the cells in the order specified in the
		// filtration. If we're in the percolation dimension + 1, re-index the
		// (d+1)-dim cells to reflect the ordering of the d-dimensional cells.
		// Otherwise, do nothing.
		if ((startDimension <= t) && (t < stopDimension)) {
			for (int j=0; j < faces; j++) column.push_back(complex->Boundary.Flat[filtration[t]][j]);
		} else {
			for (int j=0; j < faces; j++) column.push_back(remapping[complex->Boundary.Flat[t][j]]);
		}

		std::sort(column.begin(), column.end());

		boundary.set_dim(t, faces/2);
		boundary.set_col(t, column);
		column.clear();
	}

	// Compute the persistence pairs and populate a Vector to return to the
	// user. We know when the giant(/essential) cycles are born because they
	// do not appear in the list of births or deaths.
	PHATPairs pairs;
	set<int> all;
	phat::compute_persistence_pairs<PHATTwist>(pairs, boundary);
	pairs.sort();

	for (phat::index i=0; i < pairs.get_num_pairs(); i++) {
		all.insert(pairs.get_pair(i).first);
		all.insert(pairs.get_pair(i).second);
	}

	vector<int> essential;

	for (int t=0; t < filtration.size(); t++) {
		if (!all.contains(filtration[t])) essential.push_back(filtration[t]);
	}

	return essential;
}

