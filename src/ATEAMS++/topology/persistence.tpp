
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
	inline INDEX youngestFaceIndexOf(SparseVector<RingLike> cell) {
		return cell(cell.size()-1);
	}

	/**
	 * @brief 
	 */
	template <typename RingLike>
	inline void standardMatrixPolicy(
		SparseMatrix<RingLike>& Full,
		vector<int>& youngestChainLookup,
		int youngestFaceIndex,
		int markedIndex
	) {
		youngestChainLookup[youngestFaceIndex] = markedIndex;
	};

	template <typename RingLike>
	inline void twistMatrixPolicy(
		SparseMatrix<RingLike>& Full,
		vector<int>& youngestChainLookup,
		int youngestFaceIndex,
		int markedIndex
	) {
		youngestChainLookup[youngestFaceIndex] = markedIndex;
		Full.rows[youngestFaceIndex].zero();
	};

	template <typename RingLike>
	inline void reduceChains(
		SparseMatrix<RingLike>& Full,
		int start,
		int stop,
		vector<int>& youngestChainLookup,
		set<int>& marked,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options,
		function<void(SparseMatrix<RingLike>&,vector<int>&,int,int)> matrixPolicy
	) {
		SparseVector<RingLike> youngestChain;

		for (int j=start; j < stop; j++) {
			SparseVector<RingLike>& cell = Full.rows[j];

			while (cell.size() > 0 && youngestChainLookup[youngestFaceIndexOf<RingLike>(cell)] != 0) {
				youngestChain = Full.rows[youngestChainLookup[youngestFaceIndexOf<RingLike>(cell)]];
				typename RingLike::dtype q = *youngestChain.find(youngestFaceIndexOf<RingLike>(cell));

				typename RingLike::dtype s = scalar_neg(scalar_inv(q, R->ring), R->ring);
				arithmetic::SparseVectorRescaling<RingLike>(s, youngestChain, R);
				arithmetic::SparseVectorAddition<RingLike>(cell, youngestChain, R, options);
			}

			cell.compress();

			if (cell.size() > 0) {
				// Not a cell that makes the youngest chain (with which it shares
				// a face) a cycle, so we mark the index of its youngest face.
				matrixPolicy(Full, youngestChainLookup, youngestFaceIndexOf<RingLike>(cell), j);
			} else {
				// Induces a cycle, so we mark the cell.
				marked.insert(j);
			}
		}
	}

	template <typename RingLike>
	inline SparseMatrix<RingLike> reindexSparseBoundaryMatrix(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		int dimension
	) {
		// TODO this can be parallelized as well!
		// Construct an index mapping.
		vector<int> remapping(filtration.size(), 0);
		for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

		// Wow, this is going to be annoying.
		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseMatrix<RingLike> Reindexed(Full.nrow, Full.ncol);

		int startDimension = complex->Breaks[dimension][0];
		int stopDimension = complex->Breaks[dimension][1];

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

	/** @endcond */

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
	inline vector<int> standardPersistence(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		// Doing row operations on the coboundary is equivalent to column operations
		// on the boundary.
		SparseMatrix<RingLike> Full = reindexSparseBoundaryMatrix<RingLike>(complex, filtration, dimension);

		// Track which column is to be added next; track which ones are marked.
		vector<int> youngestChainSharingFace(complex->size(), 0);
		set<int> marked;

		// Top dimension of the complex; indices at which we stop and start.
		int topDimension = min(dimension+1, (int)complex->Cells.size());
		int start, stop;

		if (options.parallel->enabled) {

			// Use half the threads to work on one task, and half on the other.
			#pragma omp parallel num_threads(topDimension-dimension)
			{
				int teamSize = omp_get_max_threads()/(topDimension-dimension+1);
				omp_set_num_threads(teamSize);

				#pragma omp for
				for (int d=dimension; d <= topDimension; d++) {
					// Specify start/stop indices.
					start = complex->Breaks[d][0];
					stop = (d+1 >= complex->Cells.size()) ? complex->size() : complex->Breaks[d][1];

					options.parallel->dimensionBlocks[d].clear();

					// Reduce the chains of this dimension.
					reduceChains<RingLike>(
						Full,
						start,
						stop,
						youngestChainSharingFace,
						options.parallel->dimensionBlocks[d],
						R,
						options,
						standardMatrixPolicy<RingLike>
					);

					// options.opt->pool.detach_task(
					// 	[d, start, stop, &Full, &youngestChainSharingFace, &R, &options] {
					// 		// Flush the "marked" cells.
					// 		options.parallel->dimensionBlocks[d].clear();

					// 		// Reduce the chains of this dimension.
					// 		reduceChains<RingLike>(
					// 			Full,
					// 			start,
					// 			stop,
					// 			youngestChainSharingFace,
					// 			options.parallel->dimensionBlocks[d],
					// 			R,
					// 			options,
					// 			standardMatrixPolicy<RingLike>
					// 		);
					// 	}
					// );
				}
			}
			// #pragma omp barrier

			// options.opt->pool.wait();

			// This is constant-time, since the dimensions of the homology groups
			// are known.
			for (int d=dimension; d<= topDimension; d++) {
				for (auto k : options.parallel->dimensionBlocks[d]) {
					marked.insert(k);
				}
			}
		} else {
			for (int d=dimension; d <= topDimension; d++) {
				start = complex->Breaks[d][0];
				stop = (d+1 >= complex->Cells.size()) ? complex->size() : complex->Breaks[d][1];
				reduceChains<RingLike>(
					Full,
					start,
					stop,
					youngestChainSharingFace,
					marked,
					R,
					options,
					standardMatrixPolicy<RingLike>
				);
			}
		}

		int low = complex->Breaks[dimension][0], high = complex->Breaks[dimension][1];
		vector<int> essential;

		for (auto k : marked) {
			if (youngestChainSharingFace[k] == 0 && (low <= k && k < high)) essential.push_back(k);
		}

		return essential;
	};

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
		SparseMatrix<RingLike> Full = reindexSparseBoundaryMatrix<RingLike>(complex, filtration, dimension);

		// Track which column is to be added next; track which ones are marked.
		vector<int> youngestChainSharingFace(complex->size(), 0);
		set<int> marked;

		// Top dimension of the complex; indices at which we stop and start.
		int topDimension = min(dimension+1, (int)complex->Cells.size());
		int start, stop;

		for (int d=topDimension; d >= dimension; d--) {
			start = complex->Breaks[d][0];
			stop = (d+1 >= complex->Cells.size()) ? complex->size() : complex->Breaks[d][1];
			reduceChains<RingLike>(Full, start, stop, youngestChainSharingFace, marked, R, options, twistMatrixPolicy<RingLike>);
		}

		int low = complex->Breaks[dimension][0], high = complex->Breaks[dimension][1];
		vector<int> essential;

		for (auto k : marked) {
			if (youngestChainSharingFace[k] == 0 && (low <= k && k < high)) essential.push_back(k);
		}

		return essential;
	};

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
			if (options.parallel->enabled) essential = standardPersistence(complex, filtration, R, dimension, options);
			else essential = twistPersistence<RingLike>(complex, filtration, R, dimension, options);
		}

		std::sort(essential.begin(), essential.end());
		return essential;
	};
}

#endif
