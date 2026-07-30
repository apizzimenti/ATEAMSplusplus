
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_STANDARD_T
#define ATEAMS_TOPOLOGY_PERSISTENCE_STANDARD_T

#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_STANDARD_H
#error __FILE__ should only be included from topology/persistence/standard.h.
#endif

#include "ATEAMS++/topology/persistence/standard.h"

using namespace std;

namespace ATEAMS::topology::persistence {

	template <typename RingLike>
	inline vector<int> standard(
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

		options.serial->flush();

		// Cycle creation policy.
		auto creationPolicy = std::bind(
			policies::standardCreationPolicy<RingLike>,	// using the standard creation policy
			placeholders::_1,					// placeholder for `markedIndex`
			placeholders::_2,					// placeholder for `dim`
			std::ref(options.serial->marked)					// a reference to `marked`.
		);

		// Cycle destruction policy.
		auto destructionPolicy = std::bind(
			policies::standardDestructionPolicy<RingLike>,	// using the standard destruction policy
			placeholders::_1,						// placeholder for `cell`
			placeholders::_2,						// placeholder for `markedIndex`
			placeholders::_3,						// placeholder for `dim`,
			std::ref(options.serial->lookup)			// reference to `youngestChainLookup`.
		);

		// Reduce the blocks of the matrix.
		vector<int> endpoints = traversalPolicy(complex);

		for (int d=endpoints[0]; d <= endpoints[1]; d++) {
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				options.serial->lookup,
				d,
				R,
				policies::standardReductionPolicy<RingLike>,
				creationPolicy,
				destructionPolicy,
				options
			);
		}

		// Find essential cycles.
		return reportingPolicy(
			complex,
			options.serial->lookup,
			options.serial->marked
		);
	}


	template <typename RingLike>
	inline vector<int> standard(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		auto traversalPolicy = std::bind(
			policies::standardRestrictedTraversalPolicy<RingLike>,	// standard restricted traversal policy, since we're in a specific range
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

		return standard<RingLike>(
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
	inline vector<int> standard(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		return standard<RingLike>(
			complex,
			filtration,
			R,
			options,
			policies::fullReindexingPolicy<RingLike>,
			policies::standardFullTraversalPolicy<RingLike>,
			policies::standardFullReportingPolicy<RingLike>
		);
	}
}

#endif
