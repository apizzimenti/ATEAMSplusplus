
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
			return policies::creation::parallel<RingLike>(
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
			return policies::destruction::JIT<RingLike>(
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
			return policies::reduction::JIT<RingLike>(
				chain,
				options.parallel->lookup,
				index,
				dim,
				options.parallel->zeroed
			);
		};

		// Reduce the blocks in parallel, since they are independent of one another.
		vector<int> endpoints = traversalPolicy(complex);

		#pragma omp parallel for default(shared) schedule(static,1)
		for (int d=endpoints[0]; d >= endpoints[1]; d--) {
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
	inline vector<int> JIT(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
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
			arithmetic::ComputeOptions<RingLike>& options
		) {
			return policies::reindexing::single<RingLike>(complex, filtration, options, dimension);
		};

		return JIT<RingLike>(
			complex,
			filtration,
			R,
			options,
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
		arithmetic::ComputeOptions<RingLike>& options
	) {
		return JIT<RingLike>(
			complex,
			filtration,
			R,
			options,
			policies::reindexing::full<RingLike>,
			policies::traversal::twistFull<RingLike>,
			policies::reporting::standardFull<RingLike>
		);
	}
}

#endif
