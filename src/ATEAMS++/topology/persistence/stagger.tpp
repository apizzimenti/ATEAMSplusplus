
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
			placeholders::_1,						// placeholder for `markedIndex`
			placeholders::_2,						// placeholder for `dim`
			std::ref(options)						// a reference to `options`, for marking in parallel.
		);

		// // Cycle destruction policy.
		// auto destructionPolicy = std::bind(
		// 	policies::twistDestructionPolicy<RingLike>,	// using the twist destruction policy
		// 	placeholders::_1,					// placeholder for `cell`
		// 	placeholders::_2,					// placeholer for `markedIndex`
		// 	placeholders::_3,					// placeholder for `dim`
		// 	std::ref(options.parallel->lookup),	// reference to `youngestChainLookup`
		// 	std::ref(Full)						// reference to `Full`, for clearing
		// );

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
		
		// #pragma omp parallel default(shared)
		// {
		// 	#pragma omp for nowait schedule(static,1)
		// 	for (int d=endpoints[0]; d >= endpoints[1]; d -= 2) {
		// 		reduceBlock<RingLike>(
		// 			Full,
		// 			complex->Breaks[d],
		// 			options.parallel->lookup,
		// 			d,
		// 			R,
		// 			reductionPolicy,
		// 			creationPolicy,
		// 			destructionPolicy,
		// 			options
		// 		);
		// 	}
			
		// 	for (int d=endpoints[0]-1; d >= endpoints[1]; d -= 2) {
		// 		reduceBlock<RingLike>(
		// 			Full,
		// 			complex->Breaks[d],
		// 			options.parallel->lookup,
		// 			d,
		// 			R,
		// 			reductionPolicy,
		// 			creationPolicy,
		// 			destructionPolicy,
		// 			options
		// 		);
		// 	}
		// }

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
