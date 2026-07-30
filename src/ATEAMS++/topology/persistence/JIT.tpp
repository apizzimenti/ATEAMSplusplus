
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_JIT_T
#define ATEAMS_TOPOLOGY_PERSISTENCE_JIT_T

#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_JIT_H
#error __FILE__ should only be included from topology/persistence/JIT.h.
#endif

#include "ATEAMS++/topology/persistence/JIT.h"

using namespace std;

namespace ATEAMS::topology::persistence {
	template <typename RingLike>
	inline vector<int> JIT(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options,
		policies::ReindexingPolicy<RingLike> reindexingPolicy,
		policies::TraversalPolicy<RingLike> traversalPolicy,
		policies::ReportingPolicy<RingLike> reportingPolicy
	) {
		// Determine the endpoints and reindex the boundary matrix accordingly.
		SparseMatrix<RingLike> Full = reindexingPolicy(complex, filtration, options);

		// Flush data structures.
		options.parallel->flush();

		// Cycle creation policy.
		auto creationPolicy = std::bind(
			policies::parallelCreationPolicy<RingLike>,	// using the standard parallel creation policy
			placeholders::_1,							// placeholder for `markedIndex`
			placeholders::_2,							// placeholder for `dim`
			std::ref(options)							// a reference to `options`, for marking in parallel.
		);

		// Cycle destruction policy.
		auto destructionPolicy = std::bind(
			policies::JITDestructionPolicy<RingLike>,		// using the JIT destruction policy
			placeholders::_1,					// placeholder for `cell`
			placeholders::_2,					// placeholder for `markedIndex`
			placeholders::_3,					// placeholder for `dim`,
			std::ref(options.parallel->lookup),	// reference to `youngestChainLookup`,
			std::ref(options.parallel->zeroed)	// reference to `zeroed`.
		);

		// Reduction policy.
		auto reductionPolicy = std::bind(
			policies::JITReductionPolicy<RingLike>,		// using the JIT reduction policy
			placeholders::_1,					// placeholder for `cell`						
			placeholders::_2,					// placeholder for `lookup`
			placeholders::_3,					// placeholder for `cellIndex`
			placeholders::_4,					// placeholder for `dim`
			std::ref(options.parallel->zeroed)	// reference to `zeroed`.
		);

		// Reduce the blocks in parallel, since they are independent of one another.
		vector<int> endpoints = traversalPolicy(complex);

		#pragma omp parallel for firstprivate(Full) default(shared) schedule(static)
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
		auto traversalPolicy = std::bind(
			policies::twistRestrictedTraversalPolicy<RingLike>,	// standard restricted traversal policy, since we're in a specific range
			placeholders::_1,								// placeholder for `complex`
			dimension										// autofill the `dimension` parameter.
		);

		auto reportingPolicy = std::bind(
			policies::standardRestrictedReportingPolicy<RingLike>,	// again restricted, since we're in a range
			placeholders::_1,
			placeholders::_2,
			placeholders::_3,
			dimension
		);

		auto reindexingPolicy = std::bind(
			policies::singleReindexingPolicy<RingLike>,
			placeholders::_1,
			placeholders::_2,
			placeholders::_3,
			dimension
		);

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
			policies::fullReindexingPolicy<RingLike>,
			policies::twistFullTraversalPolicy<RingLike>,
			policies::standardFullReportingPolicy<RingLike>
		);
	}
}

#endif
