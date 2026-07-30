
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

		// // Cycle creation policy.
		// auto creationPolicy = std::bind(
		// 	policies::parallelCreationPolicy<RingLike>,	// using the standard parallel creation policy
		// 	placeholders::_1,							// placeholder for `markedIndex`
		// 	placeholders::_2,							// placeholder for `dim`
		// 	std::ref(options)							// a reference to `options`, for marking in parallel.
		// );

		auto creationPolicy = [&options](int markedIndex, int dim) {
			return policies::parallelCreationPolicy<RingLike>(markedIndex, dim, options);
		};

		// // Cycle destruction policy.
		// auto destructionPolicy = std::bind(
		// 	policies::JITDestructionPolicy<RingLike>,		// using the JIT destruction policy
		// 	placeholders::_1,					// placeholder for `cell`
		// 	placeholders::_2,					// placeholder for `markedIndex`
		// 	placeholders::_3,					// placeholder for `dim`,
		// 	std::ref(options.parallel->lookup),	// reference to `youngestChainLookup`,
		// 	std::ref(options.parallel->zeroed)	// reference to `zeroed`.
		// );

		auto destructionPolicy = [&options](SparseVector<RingLike>& chain, int markedIndex, int dim) {
			return policies::JITDestructionPolicy<RingLike>(
				chain,
				markedIndex,
				dim,
				options.parallel->lookup,
				options.parallel->zeroed
			);
		};

		// Reduction policy.
		// auto reductionPolicy = std::bind(
		// 	policies::JITReductionPolicy<RingLike>,		// using the JIT reduction policy
		// 	placeholders::_1,					// placeholder for `cell`						
		// 	placeholders::_2,					// placeholder for `lookup`
		// 	placeholders::_3,					// placeholder for `cellIndex`
		// 	placeholders::_4,					// placeholder for `dim`
		// 	std::ref(options.parallel->zeroed)	// reference to `zeroed`.
		// );

		auto reductionPolicy = [&options](SparseVector<RingLike>& chain, vector<int>& lookup, int index, int dim) {
			return policies::JITReductionPolicy<RingLike>(
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
		// auto traversalPolicy = std::bind(
		// 	policies::twistRestrictedTraversalPolicy<RingLike>,	// standard restricted traversal policy, since we're in a specific range
		// 	placeholders::_1,								// placeholder for `complex`
		// 	dimension										// autofill the `dimension` parameter.
		// );

		auto traversalPolicy = [dimension](complexes::Complex<RingLike>* complex) {
			return policies::twistRestrictedTraversalPolicy<RingLike>(complex, dimension);
		};

		// auto reportingPolicy = std::bind(
		// 	policies::standardRestrictedReportingPolicy<RingLike>,	// again restricted, since we're in a range
		// 	placeholders::_1,
		// 	placeholders::_2,
		// 	placeholders::_3,
		// 	dimension
		// );

		auto reportingPolicy = [dimension](	
			complexes::Complex<RingLike>* complex,
			vector<int>& lookup,
			set<int>& marked
		) {
			return policies::standardRestrictedReportingPolicy<RingLike>(complex, lookup, marked, dimension);
		};

		// auto reindexingPolicy = std::bind(
		// 	policies::singleReindexingPolicy<RingLike>,
		// 	placeholders::_1,
		// 	placeholders::_2,
		// 	placeholders::_3,
		// 	dimension
		// );

		auto reindexingPolicy = [dimension](
			complexes::Complex<RingLike>* complex,
			vector<int>& filtration,
			arithmetic::ComputeOptions<RingLike>& options
		) {
			return policies::singleReindexingPolicy<RingLike>(complex, filtration, options, dimension);
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
		// auto reindexingPolicy = [&](
		// 	complexes::Complex<RingLike>* complex,
		// 	vector<int>& filtration,
		// 	arithmetic::ComputeOptions<RingLike>& options
		// ) {
		// 	return helpers::reindexSparseBoundaryMatrix(complex, filtration, options);
		// };

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
