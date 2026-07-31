
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_JIT_T
#define ATEAMS_TOPOLOGY_PERSISTENCE_JIT_T

#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_JIT_H
#error __FILE__ should only be included from topology/persistence/JIT.h.
#endif

#include "ATEAMS++/topology/persistence/JIT.h"

using namespace std;

namespace ATEAMS::topology::persistence {
	template <typename RingLike>
	std::vector<int> JIT(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
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

		// Reduce the blocks in parallel, since they are independent of one another.
		vector<int> endpoints = traversalPolicy(complex);

		#pragma omp parallel for default(shared) schedule(static,1)
		for (int d=endpoints[0]; d >= endpoints[1]; d--) {
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
	inline vector<int> JIT(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		// Dimension traversal policy.
		auto traversalPolicy = [&dimension](
			complexes::Complex<RingLike>* complex
		) {
			return policies::traversal::twistRestricted<RingLike>(
				complex,
				dimension
			);
		};

		// Essential cycle reporting policy.
		auto reportingPolicy = [&dimension](	
			complexes::Complex<RingLike>* complex,
			vector<int>& lookup,
			set<int>& marked
		) {
			return policies::reporting::standardRestricted<RingLike>(
				complex,
				lookup,
				marked,
				dimension
			);
		};

		// Matrix reindexing policy.
		auto reindexingPolicy = [&dimension](
			complexes::Complex<RingLike>* complex,
			vector<int>& filtration,
			arithmetic::ComputeResources<RingLike>& resources
		) {
			return policies::reindexing::single<RingLike>(complex, filtration, resources, dimension);
		};

		return JIT<RingLike>(
			complex,
			filtration,
			R,
			resources,
			reindexingPolicy,
			traversalPolicy,
			reportingPolicy
		);
	}

	template <typename RingLike>
	inline vector<int> JIT(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		return JIT<RingLike>(
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
