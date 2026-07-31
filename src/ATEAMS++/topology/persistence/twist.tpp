
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_TWIST_T
#define ATEAMS_TOPOLOGY_PERSISTENCE_TWIST_T

#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_TWIST_H
#error __FILE__ should only be included from topology/persistence/twist.h.
#endif

#include "ATEAMS++/topology/persistence/twist.h"

using namespace std;

namespace ATEAMS::topology::persistence {
	template <typename RingLike>
	inline vector<int> twist(
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

		// Flush reusable containers.
		resources.serial->flush();

		// Cycle creation policy.
		auto creationPolicy = [&resources](
			int markedIndex,
			int dim
		) {
			policies::creation::standard<RingLike>(
				markedIndex,
				dim,
				resources.serial->marked
			);
		};


		// Cycle destruction policy.
		auto destructionPolicy = [&resources, &Full](
			SparseVector<RingLike>& chain,
			int markedIndex,
			int dim
		) {
			policies::destruction::twist<RingLike>(
				chain,
				markedIndex,
				dim,
				resources.serial->lookup,
				Full
			);
		};

		// Reduce the blocks of the matrix.
		vector<int> endpoints = traversalPolicy(complex);

		for (int d=endpoints[0]; d >= endpoints[1]; d--) {
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				resources.serial->lookup,
				d,
				R,
				policies::reduction::standard<RingLike>,
				creationPolicy,
				destructionPolicy,
				resources
			);
		}

		// Find essential cycles.
		return reportingPolicy(
			complex,
			resources.serial->lookup,
			resources.serial->marked
		);
	};


	template <typename RingLike>
	inline vector<int> twist(
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

		return twist<RingLike>(
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
	inline vector<int> twist(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		return twist<RingLike>(
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
