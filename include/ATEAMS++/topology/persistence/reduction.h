
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

	/**
	 * @brief Template specialization for reducing Z/2Z rows. This won't get used,
	 * 	but we have to have it here for everything to build properly.
	 */
	template <>
	inline void reduceChain<Z2>(
		SparseMatrix<Z2>& Full,
		std::vector<int>& lookup,
		SparseVector<Z2>& cell,
		Ring* R,
		arithmetic::ComputeResources<Z2>& resources
	) {
		// Get the youngest chain (column) sharing a face with `cell` and zero
		// out the pivot row.
		SparseVector<Z2> youngest = Full.rows[lookup[helpers::youngestOf<Z2>(cell)]];
		arithmetic::SparseVectorAddition<Z2>(cell, youngest, R, resources);
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

