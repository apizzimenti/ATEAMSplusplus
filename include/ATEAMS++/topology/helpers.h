
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_HELPERS_H
#define AETAMS_TOPOLOGY_PERSISTENCE_HELPERS_H

#include <set>
#include <vector>

/** @cond */
namespace ATEAMS::topology::helpers {

	template <typename RingLike>
	inline INDEX youngestOf(SparseVector<RingLike>& cell) {
		return cell(cell.size()-1);
	}

	template <typename RingLike>
	inline SparseMatrix<RingLike> reindexSparseBoundaryMatrix(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		arithmetic::ComputeResources<RingLike>& options
	) {
		// Construct an index mapping.
		vector<int> remapping(filtration.size());
		
		#pragma omp parallel for shared(filtration, remapping) if(options.parallel->enabled)
		for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

		// Duplicate the matrix.
		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseMatrix<RingLike> Reindexed(Full.nrow, Full.ncol);

		#pragma omp parallel for shared(filtration, remapping)
		for (int t=0; t < Full.nrow; t++) {
			SparseVector<RingLike>& row = Reindexed.rows[t];
			SparseVector<RingLike> orow = Full.rows[filtration[t]];

			for (int i=0; i < orow.size(); i++) {
				row.push_back(
					(INDEX)remapping[orow(i)],
					(typename RingLike::dtype)orow[i]
				);
			}
		}

		Reindexed.compress();
		return Reindexed;
	}

	
	template <typename RingLike>
	inline SparseMatrix<RingLike> reindexSparseBoundaryMatrix(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		int dimension,
		arithmetic::ComputeResources<RingLike>& options
	) {
		// Construct an index mapping.
		vector<int> remapping(filtration.size());
		
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
}
/** @endcond */

#endif