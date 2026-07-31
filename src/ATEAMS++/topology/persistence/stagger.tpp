
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_STAGGER_T
#define ATEAMS_TOPOLOGY_PERSISTENCE_STAGGER_T

#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_STAGGER_H
#error __FILE__ should only be included from topology/persistence/stagger.h.
#endif

#include "ATEAMS++/topology/persistence/stagger.h"
#include "ATEAMS++/topology/persistence/twist.h"

using namespace std;

namespace ATEAMS::topology::persistence {
	template <typename RingLike>
	inline vector<int> stagger(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources,
		auto& reindexingPolicy,
		auto& traversalPolicy,
		auto& reportingPolicy
	) {
		// Determine the endpoints and reindex the boundary matrix accordingly.
		SparseMatrix<RingLike> Full = reindexingPolicy(complex, filtration, resources);

		// Flush data structures.
		resources.parallel->flush();

		// Cycle creation policy.
		auto creationPolicy = [&resources](
			int markedIndex,
			int dim
		) {
			return policies::creation::parallel<RingLike>(
				markedIndex,
				dim,
				resources
			);
		};

		// Cycle destruction policy.
		auto destructionPolicy = [&resources](
			SparseVector<RingLike>& chain,
			int markedIndex,
			int dim
		) {
			return policies::destruction::JIT<RingLike>(
				chain,
				markedIndex,
				dim,
				resources.parallel->lookup,
				resources.parallel->zeroed
			);
		};

		// Reduction policy.
		auto reductionPolicy = [&resources](
			SparseVector<RingLike>& chain,
			vector<int>& lookup,
			int index,
			int dim
		) {
			return policies::reduction::JIT<RingLike>(
				chain,
				resources.parallel->lookup,
				index,
				dim,
				resources.parallel->zeroed
			);
		};

		// Stagger the block reductions, so we can incorporate some of the
		// clearing optimization benefits.
		vector<int> endpoints = traversalPolicy(complex);

		#pragma omp parallel for default(shared) schedule(static,1)
		for (int d=endpoints[0]; d >= endpoints[1]; d -= 2) {
			// Reduce the dth block,
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				resources.parallel->lookup,
				d,
				R,
				reductionPolicy,
				creationPolicy,
				destructionPolicy,
				resources
			);

			// Then the d-1th block, if d-1 >= 0.
			if (d-1 >= 0) {
				reduceBlock<RingLike>(
					Full,
					complex->Breaks[d-1],
					resources.parallel->lookup,
					d-1,
					R,
					reductionPolicy,
					creationPolicy,
					destructionPolicy,
					resources
				);
			}
		}

		// Re-constitute the marked columns.
		set<int> marked;
		for (int d=endpoints[0]; d >= endpoints[1]; d--) {
			for (auto& k : resources.parallel->marked[d]) marked.insert(k);
		}

		// Find essential cycles.
		return reportingPolicy(
			complex,
			resources.parallel->lookup,
			marked
		);
	}


	template <typename RingLike>
	inline vector<int> stagger(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		return twist<RingLike>(
			complex,
			filtration,
			R,
			dimension,
			resources
		);
	}

	template <typename RingLike>
	inline vector<int> stagger(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		return stagger<RingLike>(
			complex,
			filtration,
			R,
			resources,
			policies::reindexing::full<RingLike>,
			policies::traversal::twistFull<RingLike>,
			policies::reporting::standardFull<RingLike>
		);
	}
}

#endif
