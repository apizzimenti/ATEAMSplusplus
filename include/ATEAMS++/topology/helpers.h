
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
		arithmetic::ComputeResources<RingLike>& resources
	) {
		// Construct an index mapping.
		if (resources.parallel->remapping.size() != filtration.size()) {
			resources.parallel->remapping = vector<int>(filtration.size(), 0);
		}

		#pragma omp parallel for default(shared) if(resources.parallel->enabled && filtration.size() > 1e5)
		for (int t=0; t < filtration.size(); t++) resources.parallel->remapping[filtration[t]] = t;

		// resources.opt->pool.detach_loop(
		// 	0, filtration.size(),
		// 	[&](const size_t t) {
		// 		resources.parallel->remapping[filtration[t]] = t;
		// 	}
		// );
		// resources.opt->pool.wait();

		// Duplicate the matrix.
		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseMatrix<RingLike> Reindexed(Full.nrow, Full.ncol);

		#pragma omp parallel for default(shared) if(resources.parallel->enabled && filtration.size() > 1e5)
		for (int t=0; t < Full.nrow; t++) {
			SparseVector<RingLike>& row = Reindexed.rows[t];
			SparseVector<RingLike> orow = Full.rows[filtration[t]];

			for (int i=0; i < orow.size(); i++) {
				row.push_back(
					(INDEX)resources.parallel->remapping[orow(i)],
					(typename RingLike::dtype)orow[i]
				);
			}
		}

		// If we haven't yet decided on a preferred number of blocks, let's
		// compute that now.
		// if (resources.typed->blocks < 0) {
		// 	resources.typed->defaultPreferredBlocks(Full.nrow);
		// 	resources.typed->build();
		// }

		// resources.opt->pool.detach_blocks(
		// 	0, Full.nrow,
		// 	[&](const size_t start, const size_t end) {
		// 		for (int t=start; t < end; t++) {
		// 			SparseVector<RingLike>& row = Reindexed.rows[t];
		// 			SparseVector<RingLike> orow = Full.rows[filtration[t]];

		// 			for (int i=0; i < orow.size(); i++) {
		// 				row.push_back(
		// 					(INDEX)resources.parallel->remapping[orow(i)],
		// 					(typename RingLike::dtype)orow[i]
		// 				);
		// 			}
		// 		}
		// 	},
		// 	resources.typed->blocks
		// );
		// resources.opt->pool.wait();

		Reindexed.compress();
		return Reindexed;
	}

	
	template <typename RingLike>
	inline SparseMatrix<RingLike> reindexSparseBoundaryMatrix(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		int dimension,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		// Construct an index mapping.
		if (resources.parallel->remapping.size() != filtration.size()) {
			resources.parallel->remapping = vector<int>(filtration.size(), 0);
		}
		
		#pragma omp parallel for default(shared) if(resources.parallel->enabled && filtration.size() > 1e5)
		for (int t=0; t < filtration.size(); t++) resources.parallel->remapping[filtration[t]] = t;

		// resources.opt->pool.detach_loop(
		// 	0, filtration.size(),
		// 	[&](const size_t t) {
		// 		resources.parallel->remapping[filtration[t]] = t;
		// 	}
		// );
		// resources.opt->pool.wait();

		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseMatrix<RingLike> Reindexed(Full.nrow, Full.ncol);

		int startDimension = complex->Breaks[dimension][0];
		int stopDimension = complex->Breaks[dimension][1];

		#pragma omp parallel for default(shared) if (resources.parallel->enabled && filtration.size() > 1e5)
		for (int t=0; t < Full.nrow; t++) {
			if ((startDimension <= t) && (t < stopDimension)) {
				Reindexed.rows[t] = Full.rows[filtration[t]];
			} else {
				SparseVector<RingLike>& row = Reindexed.rows[t];
				SparseVector<RingLike> orow = Full.rows[t];

				for (int i=0; i < orow.size(); i++) {
					row.push_back(
						(INDEX)resources.parallel->remapping[orow(i)],
						(typename RingLike::dtype)orow[i]
					);
				}
			}
		}

		// If we haven't yet decided on a preferred number of blocks, let's
		// compute that now.
		// if (resources.typed->blocks < 0) {
		// 	resources.typed->defaultPreferredBlocks(Full.nrow);
		// 	resources.typed->build();
		// }

		// resources.opt->pool.detach_blocks(
		// 	0, Full.nrow,
		// 	[&](const size_t start, const size_t end) {
		// 		if ((startDimension <= t) && (t < stopDimension)) {
		// 			Reindexed.rows[t] = Full.rows[filtration[t]];
		// 		} else {
		// 			SparseVector<RingLike>& row = Reindexed.rows[t];
		// 			SparseVector<RingLike> orow = Full.rows[t];

		// 			for (int i=0; i < orow.size(); i++) {
		// 				row.push_back(
		// 					(INDEX)resources.parallel->remapping[orow(i)],
		// 					(typename RingLike::dtype)orow[i]
		// 				);
		// 			}
		// 		}
		// 	},
		// 	resources.typed->blocks
		// );
		// resources.opt->pool.wait();

		Reindexed.compress();
		return Reindexed;
	}
}
/** @endcond */

#endif