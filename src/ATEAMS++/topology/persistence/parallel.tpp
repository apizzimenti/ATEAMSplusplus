
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_PARALLEL_T
#define ATEAMS_TOPOLOGY_PERSISTENCE_PARALLEL_T

#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_PARALLEL_H
#error __FILE__ should only be included from topology/persistence/parallel.h.
#endif

#include "ATEAMS++/topology/persistence/parallel.h"

using namespace std;

namespace ATEAMS::topology::persistence {
	template <typename RingLike>
	inline vector<int> parallel(
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
			return policies::destruction::standard<RingLike>(
				chain,
				markedIndex,
				dim,
				options.parallel->lookup
			);
		};

		// Reduce the blocks in parallel, since they are independent of one another.
		vector<int> endpoints = traversalPolicy(complex);

		#pragma omp parallel for firstprivate(Full) default(shared) schedule(static,1)
		for (int d=endpoints[0]; d <= endpoints[1]; d++) {
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				options.parallel->lookup,
				d,
				R,
				policies::reduction::standard<RingLike>,
				creationPolicy,
				destructionPolicy,
				options
			);
		}

		// Re-constitute the marked columns.
		set<int> marked;
		for (int d=endpoints[0]; d <= endpoints[1]; d++) {
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
	inline vector<int> parallel(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		auto traversalPolicy = [&dimension](
			complexes::Complex<RingLike>* complex
		) {
			return policies::traversal::standardRestricted<RingLike>(
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

		return parallel<RingLike>(
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
	inline vector<int> parallel(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		return parallel<RingLike>(
			complex,
			filtration,
			R,
			options,
			policies::reindexing::full<RingLike>,
			policies::traversal::standardFull<RingLike>,
			policies::reporting::standardFull<RingLike>
		);
	}
}

#endif
