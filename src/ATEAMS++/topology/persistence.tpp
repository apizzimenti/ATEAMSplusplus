
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_T
#define ATEAMS_TOPOLOGY_PERSISTENCE_T

#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_H
#error __FILE__ should only be included from topology/persistence.h.
#endif

#include "ATEAMS++/topology/persistence.h"
#include "ATEAMS++/arithmetic/arithmetic.h"

#include <phat/compute_persistence_pairs.h>
#include <phat/helpers/misc.h>

using namespace SparseRREF;
using namespace std;

namespace ATEAMS::topology {
	/** @cond */
	template <typename RingLike>
	inline INDEX youngestOf(SparseVector<RingLike> cell) {
		return cell(cell.size()-1);
	}

	template <typename RingLike>
	inline SparseMatrix<RingLike> reindexSparseBoundaryMatrix(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		// Construct an index mapping.
		vector<int> remapping(filtration.size(), 0);
		
		#pragma omp parallel for shared(filtration, remapping) if(options.parallel->enabled)
		for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseMatrix<RingLike> Reindexed(Full.nrow, Full.ncol);

		int startDimension = complex->Breaks[dimension][0];
		int stopDimension = complex->Breaks[dimension][1];

		#pragma omp parallel for shared(Full, Reindexed) if (options.parallel->enabled)
		for (int t=0; t < Full.nrow; t++) {
			if ((startDimension <= t) && (t < stopDimension)) {
				Reindexed.rows[t] = Full.rows[filtration[t]];
			} else {
				SparseVector<RingLike>& row = Reindexed.rows[t];
				SparseVector<RingLike> orow = Full.rows[t];

				for (int i=0; i < orow.size(); i++) {
					row.push_back(
						(INDEX)remapping[orow(i)],
						(typename RingLike::dtype)orow[i]
					);
				}
			}
		}

		Reindexed.compress();
		return Reindexed;
	}


	template <typename RingLike>
	inline bool standardReductionPolicy(
		SparseVector<RingLike>& cell,
		vector<int>& lookup,
		int dim
	) {
		return (cell.size() > 0) && (lookup[youngestOf<RingLike>(cell)] != 0);
	}


	template <typename RingLike>
	inline void standardCreationPolicy(
		int markedIndex,
		int dim,
		set<int>& marked
	) {
		marked.insert(markedIndex);
	}


	template <typename RingLike>
	inline void standardDestructionPolicy(
		SparseVector<RingLike>& cell,
		int markedIndex,
		int dim,
		vector<int>& lookup
	) {
		lookup[youngestOf<RingLike>(cell)] = markedIndex;
	}


	inline vector<int> standardReportingPolicy(
		int low,
		int high,
		vector<int>& lookup,
		set<int>& marked
	) {
		vector<int> essential;

		for (int k : marked) {
			if (lookup[k] == 0 && (low <= k && k < high)) essential.push_back(k);
		}

		return essential;
	}


	template <typename RingLike>
	inline void twistDestructionPolicy(
		SparseVector<RingLike>& cell,
		int markedIndex,
		int dim,
		vector<int>& lookup,
		SparseMatrix<RingLike>& Full
	) {
		lookup[youngestOf<RingLike>(cell)] = markedIndex;
		Full.rows[youngestOf<RingLike>(cell)].zero();
	};


	template <typename RingLike>
	inline void standardParallelCreationPolicy(
		int markedIndex,
		int dim,
		set<int>& marked,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		options.parallel->markedByThread[dim].insert(markedIndex);
	}


	template <typename RingLike>
	inline void reduceChain(
		SparseMatrix<RingLike>& Full,
		vector<int>& lookup,
		SparseVector<RingLike>& cell,
		Ring* R
	) {
		// Get the youngest chain (column) sharing a face with `cell` and zero
		// out the pivot row.
		SparseVector<RingLike> youngest = Full.rows[lookup[youngestOf<RingLike>(cell)]];

		typename RingLike::dtype q = *youngest.find(youngestOf<RingLike>(cell));
		typename RingLike::dtype s = scalar_neg(scalar_inv(q, R->ring), R->ring);

		arithmetic::SparseVectorRescaling<RingLike>(s, youngest, R);
		arithmetic::SparseVectorAddition<RingLike>(cell, youngest, R);
	}


	template <typename RingLike>
	inline void reduceBlock(
		SparseMatrix<RingLike>& Full,
		vector<int>& endpoints,
		vector<int>& lookup,
		int dim,
		Ring *R,
		function<bool(SparseVector<RingLike>&,vector<int>&,int)> reductionPolicy,
		function<void(int,int)> creationPolicy,
		function<void(SparseVector<RingLike>&,int,int)> destructionPolicy,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		for (int j=endpoints[0]; j < endpoints[1]; j++) {
			SparseVector<RingLike>& cell = Full.rows[j];

			while (reductionPolicy(cell, lookup, dim)) reduceChain<RingLike>(Full, lookup, cell, R);

			if (cell.size() > 0) destructionPolicy(cell, j, dim);
			else creationPolicy(j, dim);
		}
	}

	/** @endcond */


	template <typename RingLike>
	inline vector<int> standardPersistence(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		SparseMatrix<RingLike> Full = reindexSparseBoundaryMatrix<RingLike>(complex, filtration, dimension, options);

		// Track which column is to be added next; track which ones are marked.
		vector<int> youngestChainLookup(complex->size(), 0);
		set<int> marked;

		// Top dimension of the complex; indices at which we stop and start.
		int topDimension = min(dimension+1, (int)complex->Cells.size());

		// Cycle creation policy.
		auto creationPolicy = std::bind(
			standardCreationPolicy<RingLike>,	// using the standard creation policy
			placeholders::_1,					// placeholder for `markedIndex`
			placeholders::_2,					// placeholder for `dim`
			std::ref(marked)					// a reference to `marked`.
		);

		// Cycle destruction policy.
		auto destructionPolicy = std::bind(
			standardDestructionPolicy<RingLike>,	// using the standard destruction policy
			placeholders::_1,						// placeholder for `cell`
			placeholders::_2,						// placeholder for `markedIndex`
			placeholders::_3,						// placeholder for `dim`,
			std::ref(youngestChainLookup)			// reference to `youngestChainLookup`.
		);


		// Reduce the blocks of the matrix.
		for (int d=dimension; d <= topDimension; d++) {
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				youngestChainLookup,
				d,
				R,
				standardReductionPolicy<RingLike>,
				creationPolicy,
				destructionPolicy,
				options
			);
		}

		// Find essential cycles.
		return standardReportingPolicy(
			complex->Breaks[dimension][0],
			complex->Breaks[dimension][1],
			youngestChainLookup,
			marked
		);
	}


	template <typename RingLike>
	inline vector<int> twistPersistence(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		// Doing row operations on the coboundary is equivalent to column operations
		// on the boundary.
		SparseMatrix<RingLike> Full = reindexSparseBoundaryMatrix<RingLike>(complex, filtration, dimension, options);

		// Track which column is to be added next; track which ones are marked.
		vector<int> youngestChainLookup(complex->size(), 0);
		set<int> marked;

		// Top dimension of the complex; indices at which we stop and start.
		int topDimension = min(dimension+1, (int)complex->Cells.size());

		// Cycle creation policy.
		auto creationPolicy = std::bind(
			standardCreationPolicy<RingLike>,	// using the standard creation policy
			placeholders::_1,					// placeholder for `markedIndex`
			placeholders::_2,					// placeholder for `dim`
			std::ref(marked)					// a reference to `marked`.
		);

		// Cycle destruction policy.
		auto destructionPolicy = std::bind(
			twistDestructionPolicy<RingLike>,	// using the twist destruction policy
			placeholders::_1,					// placeholder for `cell`
			placeholders::_2,					// placeholer for `markedIndex`
			placeholders::_3,					// placeholder for `dim`
			std::ref(youngestChainLookup),		// reference to `youngestChainLookup`
			std::ref(Full)						// reference to `Full`, for clearing
		);

		// Reduce the blocks of the matrix.
		for (int d=topDimension; d >= dimension; d--) {
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				youngestChainLookup,
				d,
				R,
				standardReductionPolicy<RingLike>,
				creationPolicy,
				destructionPolicy,
				options
			);
		}

		// Find essential cycles.
		return standardReportingPolicy(
			complex->Breaks[dimension][0],
			complex->Breaks[dimension][1],
			youngestChainLookup,
			marked
		);
	};


	template <typename RingLike>
	inline vector<int> standardParallelPersistence(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		SparseMatrix<RingLike> Full = reindexSparseBoundaryMatrix<RingLike>(complex, filtration, dimension, options);

		// Track which column is to be added next; track which ones are marked.
		vector<int> youngestChainLookup(complex->size(), 0);
		set<int> marked;

		// Top dimension of the complex; indices at which we stop and start.
		int topDimension = min(dimension+1, (int)complex->Cells.size());

		// Cycle creation policy.
		auto creationPolicy = std::bind(
			standardParallelCreationPolicy<RingLike>,	// using the standard parallel creation policy
			placeholders::_1,							// placeholder for `markedIndex`
			placeholders::_2,							// placeholder for `dim`
			std::ref(marked),							// a dummy reference to `marked`
			std::ref(options)							// a reference to `options`, for marking in parallel.
		);

		// Cycle destruction policy.
		auto destructionPolicy = std::bind(
			standardDestructionPolicy<RingLike>,	// using the standard destruction policy
			placeholders::_1,						// placeholder for `cell`
			placeholders::_2,						// placeholder for `markedIndex`
			placeholders::_3,						// placeholder for `dim`,
			std::ref(youngestChainLookup)			// reference to `youngestChainLookup`.
		);

		// Reduce the blocks in parallel, since they are independent of one another.
		#pragma omp parallel for firstprivate(Full) shared(youngestChainLookup)
		for (int d=dimension; d <= topDimension; d++) {
			options.parallel->flush(d);

			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				youngestChainLookup,
				d,
				R,
				standardReductionPolicy<RingLike>,
				creationPolicy,
				destructionPolicy,
				options
			);
		}

		// Re-constitute the marked columns.
		for (int d=dimension; d <= topDimension; d++) {
			for (auto& k : options.parallel->markedByThread[d]) marked.insert(k);
		}

		// Find essential cycles.
		return standardReportingPolicy(
			complex->Breaks[dimension][0],
			complex->Breaks[dimension][1],
			youngestChainLookup,
			marked
		);
	}

	template <typename RingLike>
	inline vector<int> PHATPersistence(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		int dimension
	) {
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

		// TODO optimize, this takes way too long
		vector<int> remapping(filtration.size(), 0);
		for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

		PHATBoundaryMatrix boundary;
		PHATColumn column;
		int faces;

		int startDimension = complex->Breaks[dimension][0];
		int stopDimension = complex->Breaks[dimension][1];

		boundary.set_num_cols(complex->Boundary.Flat.size());

		// this also takes an insane amount of time; no need to go all the way
		// up.
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

		int low = complex->Breaks[dimension][0], high = complex->Breaks[dimension][1];
		vector<int> essential;

		for (int t=0; t < filtration.size(); t++) {
			if (!all.contains(filtration[t]) && (low <= filtration[t] && filtration[t] < high)){
				essential.push_back(filtration[t]);
			}
		}

		return essential;
	}


	template <typename RingLike>
	std::vector<int> persistence(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		vector<int> essential;

		if (R->characteristic < 3) essential = PHATPersistence<RingLike>(complex, filtration, dimension);
		else {
			if (options.parallel->enabled) essential = standardParallelPersistence(complex, filtration, R, dimension, options);
			else essential = twistPersistence<RingLike>(complex, filtration, R, dimension, options);
		}

		std::sort(essential.begin(), essential.end());
		return essential;
	};
}

#endif
