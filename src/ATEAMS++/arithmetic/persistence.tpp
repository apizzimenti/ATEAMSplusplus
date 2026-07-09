
#ifndef ATEAMS_ARITHMETIC_PERSISTENCE_T
#define ATEAMS_ARITHMETIC_PERSISTENCE_T

#ifndef ATEAMS_ARITHMETIC_PERSISTENCE_H
#error __FILE__ should only be included from arithmetic/persistence.h.
#endif

#include "ATEAMS++/common.h"
#include "ATEAMS++/util.h"
#include "ATEAMS++/complexes/Complex.h"
#include "ATEAMS++/arithmetic/persistence.h"
#include "ATEAMS++/arithmetic/options.h"

#include <phat/compute_persistence_pairs.h>

#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_vec.h>
#include <SparseRREF/scalar.h>

using namespace SparseRREF;
using namespace std;

namespace ATEAMS {
	namespace arithmetic {
		/** @cond */

		// Type specialization for adding sparse vectors.
		void add_sparse_vectors(SparseVector<Zp> &A, SparseVector<Zp> B, Zp R) {
			sparse_vec_add<index_t>(A, B, R.ring);
		}

		void add_sparse_vectors(SparseVector<Q> &A, SparseVector<Q> B, Q R) {
			sfmpq_vec_addsub_mul<index_t,false>(A, B, (Q::dtype)1);
		}


		template <typename RingLike>
		index_t youngestFaceIndexOf(SparseVector<RingLike> cell) {
			return cell(cell.size()-1);
		}


		template <typename RingLike>
		SparseMatrix<RingLike> arithmetic::reindexSparseBoundaryMatrix(
			ATEAMS::complexes::Complex<RingLike>* complex, vector<int> filtration, int dimension
		) {
			// Construct an index mapping.
			map<int,int> remapping;
			for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

			// Wow, this is going to be annoying.
			SparseMatrix<RingLike> Full = complex->Coboundary.Full;
			SparseMatrix<RingLike> Reindexed(Full.nrow, Full.ncol);

			int startDimension = complex->breaks[dimension][0];
			int stopDimension = complex->breaks[dimension][1];

			for (int t=0; t < Full.nrow; t++) {
				if ((startDimension <= t) && (t < stopDimension)) {
					Reindexed.rows[t] = Full.rows[filtration[t]];
				} else {
					SparseVector<RingLike>& row = Reindexed.rows[t];
					SparseVector<RingLike> orow = Full.rows[t];

					for (int i=0; i < orow.size(); i++) {
						row.push_back(
							(index_t)remapping[orow(i)],
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
		vector<int> arithmetic::PHATPersistence(
			ATEAMS::complexes::Complex<RingLike>* complex,
			vector<int> filtration,
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


		template <typename RingLike>
		vector<int> arithmetic::TwistPersistence(
			ATEAMS::complexes::Complex<RingLike>* complex,
			vector<int> filtration,
			Ring* R,
			int dimension
		) {
			// Doing row operations on the coboundary is equivalent to column operations
			// on the boundary.
			SparseMatrix<RingLike> Full = reindexSparseBoundaryMatrix<RingLike>(complex, filtration, dimension);
			SparseVector<RingLike> youngestFace;

			// Track which column is to be added next; track which ones are marked.
			vector<int> nextColumnAdded(complex->size(), 0);
			set<int> marked;

			// Top dimension of the complex; indices at which we stop and start.
			int topDimension = complex->Cells.size(), lowestCellIndex, highestCellIndex;
			typename RingLike::dtype youngestFaceCoefficient;

			for (int d=topDimension-1; d > dimension-1; d--) {
				lowestCellIndex = complex->breaks[d][0];
				highestCellIndex = (d+1 >= complex->Cells.size()) ? complex->size() : complex->breaks[d][1];

				for (int j=lowestCellIndex; j < highestCellIndex; j++) {
					SparseVector<RingLike>& cell = Full.rows[j];

					while (cell.size() > 0 && nextColumnAdded[youngestFaceIndexOf<RingLike>(cell)] != 0) {
						youngestFace = Full.rows[nextColumnAdded[youngestFaceIndexOf<RingLike>(cell)]];
						youngestFaceCoefficient = *youngestFace.find(youngestFaceIndexOf<RingLike>(cell));

						// TODO parallelization stuff here? Can't go across columns, so maybe
						// within the column? SparseRREF/FLINT probably do that already tho.
						sparse_vec_rescale<index_t,typename RingLike::dtype>(
							youngestFace,
							scalar_neg(scalar_inv(youngestFaceCoefficient, R->ring), R->ring),
							R->ring
						);
						add_sparse_vectors(cell, youngestFace, R->ring);

						cell.compress();
					}

					if (cell.size() > 0) {
						nextColumnAdded[youngestFaceIndexOf<RingLike>(cell)] = j;
						Full.rows[youngestFaceIndexOf<RingLike>(cell)].zero();
					} else {
						marked.insert(j);
					}
				}
			}

			int low = complex->breaks[dimension][0], high = complex->breaks[dimension][1];
			vector<int> essential;

			for (auto k : marked) {
				if (nextColumnAdded[k] == 0 && (low <= k && k < high)) essential.push_back(k);
			}

			return essential;
		}
	}
}

#endif
