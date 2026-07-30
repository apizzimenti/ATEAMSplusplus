
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
		auto& reindexingPolicy,
		auto& traversalPolicy,
		auto& reportingPolicy
	) {
		// Determine the endpoints and reindex the boundary matrix accordingly.
		SparseMatrix<RingLike> Full = reindexingPolicy(complex, filtration, options);

		// Flush reusable serial containers.
		options.serial->flush();

		// Cycle creation policy.
		auto creationPolicy = [&options](
			int markedIndex,
			int dim
		) {
			policies::standardCreationPolicy<RingLike>(
				markedIndex,
				dim,
				options.serial->marked
			);
		};

		// Cycle destruction policy.
		auto destructionPolicy = [&options](
			SparseVector<RingLike>& chain,
			int markedIndex,
			int dim
		) {
			policies::standardDestructionPolicy<RingLike>(
				chain,
				markedIndex,
				dim,
				options.serial->lookup
			);
		};

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
		auto traversalPolicy = [&dimension](
			complexes::Complex<RingLike>* complex
		) {
			return policies::standardRestrictedTraversalPolicy<RingLike>(
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
			return policies::standardRestrictedReportingPolicy<RingLike>(
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
			return policies::singleReindexingPolicy<RingLike>(complex, filtration, options, dimension);
		};

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
