
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_REDUCTION_H
#define ATEAMS_TOPOLOGY_PERSISTENCE_REDUCTION_H

#include "ATEAMS++/topology/persistence/policies.h"

#include <vector>
#include <set>


/** @brief Persistent homology routines. */
namespace ATEAMS::topology::persistence {
	/** @cond */
	template <typename RingLike>
	inline void reduceChain(
		SparseMatrix<RingLike>& Full,
		std::vector<int>& lookup,
		SparseVector<RingLike>& cell,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		// Get the youngest chain (column) sharing a face with `cell` and zero
		// out the pivot row.
		SparseVector<RingLike> youngest = Full.rows[lookup[helpers::youngestOf<RingLike>(cell)]];

		typename RingLike::dtype q = *youngest.find(helpers::youngestOf<RingLike>(cell));
		typename RingLike::dtype s = resources.arithmetic->negate[resources.arithmetic->invert[q]];

		arithmetic::SparseVectorRescaling<RingLike>(s, youngest, R, resources);
		arithmetic::SparseVectorAddition<RingLike>(cell, youngest, R, resources);
	}


	template <typename RingLike>
	inline void reduceBlock(
		SparseMatrix<RingLike>& Full,
		std::vector<int>& endpoints,
		std::vector<int>& lookup,
		int dim,
		Ring *R,
		auto& reductionPolicy,
		auto& creationPolicy,
		auto& destructionPolicy,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		for (int j=endpoints[0]; j < endpoints[1]; j++) {
			SparseVector<RingLike>& cell = Full.rows[j];

			while (reductionPolicy(cell, lookup, j, dim)) reduceChain<RingLike>(Full, lookup, cell, R, resources);

			if (cell.size() > 0) destructionPolicy(cell, j, dim);
			else creationPolicy(j, dim);
		}
	}
	/** @endcond */
}

#endif

