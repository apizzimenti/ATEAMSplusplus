
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
		arithmetic::ComputeOptions<RingLike>& options,
		auto& reindexingPolicy,
		auto& traversalPolicy,
		auto& reportingPolicy
	) {
		// Determine the endpoints and reindex the boundary matrix accordingly.
		SparseMatrix<RingLike> Full = reindexingPolicy(complex, filtration, options);

		// Flush data structures.
		options.parallel->flush();

		// Cycle creation policy.
		auto creationPolicy = [&options](
			int markedIndex,
			int dim
		) {
			return policies::parallelCreationPolicy<RingLike>(
				markedIndex,
				dim,
				options
			);
		};

		// Cycle destruction policy.
		auto destructionPolicy = [&options](
			SparseVector<RingLike>& chain,
			int markedIndex,
			int dim
		) {
			return policies::JITDestructionPolicy<RingLike>(
				chain,
				markedIndex,
				dim,
				options.parallel->lookup,
				options.parallel->zeroed
			);
		};

		// Reduction policy.
		auto reductionPolicy = [&options](
			SparseVector<RingLike>& chain,
			vector<int>& lookup,
			int index,
			int dim
		) {
			return policies::JITReductionPolicy<RingLike>(
				chain,
				options.parallel->lookup,
				index,
				dim,
				options.parallel->zeroed
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
				options.parallel->lookup,
				d,
				R,
				reductionPolicy,
				creationPolicy,
				destructionPolicy,
				options
			);

			// Then the d-1th block, if d-1 >= 0.
			if (d-1 >= 0) {
				reduceBlock<RingLike>(
					Full,
					complex->Breaks[d-1],
					options.parallel->lookup,
					d-1,
					R,
					reductionPolicy,
					creationPolicy,
					destructionPolicy,
					options
				);
			}
		}

		// Re-constitute the marked columns.
		set<int> marked;
		for (int d=endpoints[0]; d >= endpoints[1]; d--) {
			for (auto& k : options.parallel->marked[d]) marked.insert(k);
		}

		// Find essential cycles.
		return reportingPolicy(
			complex,
			options.parallel->lookup,
			marked
		);
	}


	template <typename RingLike>
	inline vector<int> stagger(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		return twist<RingLike>(
			complex,
			filtration,
			R,
			dimension,
			options
		);
	}

	template <typename RingLike>
	inline vector<int> stagger(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		return stagger<RingLike>(
			complex,
			filtration,
			R,
			options,
			policies::fullReindexingPolicy<RingLike>,
			policies::twistFullTraversalPolicy<RingLike>,
			policies::standardFullReportingPolicy<RingLike>
		);
	}
}

#endif
