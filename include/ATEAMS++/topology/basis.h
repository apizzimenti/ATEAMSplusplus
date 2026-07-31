
#ifndef ATEAMS_TOPOLOGY_BASIS_H
#define ATEAMS_TOPOLOGY_BASIS_H


#include "ATEAMS++/topology/persistence/policies.h"
#include "ATEAMS++/topology/persistence/reduction.h"


namespace ATEAMS::topology {
	template <typename RingLike>
	inline SparseBases<RingLike> basis(
		complexes::Complex<RingLike>* complex,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		// Get the full boundary matrix; no need to reindex.
		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseBases<RingLike> bases(complex->Cells.size(), SparseBasis<RingLike>());

		// Flush reusable containers.
		resources.serial->flush();

		// Cycle creation policy.
		auto creationPolicy = [&resources](
			int index,
			int dim
		) {
			persistence::policies::creation::standard<RingLike>(
				index,
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
			persistence::policies::destruction::twist<RingLike>(
				chain,
				markedIndex,
				dim,
				resources.serial->lookup,
				Full
			);
		};

		// Reduce blocks.
		vector<int> endpoints = persistence::policies::traversal::twistFull(complex);

		for (int d=endpoints[0]; d >= endpoints[1]; d--) {
			persistence::reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				resources.serial->lookup,
				d,
				R,
				persistence::policies::reduction::standard<RingLike>,
				creationPolicy,
				destructionPolicy,
				resources
			);
		}



		// Find essential cycles.
		for (int k : resources.serial->marked) {
			if (resources.serial->lookup[k] == 0) {
				int dim;

				// Determine the dimension of the cell.
				for (int d=0; d < complex->Breaks.size(); d++) {
					int lo = complex->Breaks[d][0];
					int hi = complex->Breaks[d][1];

					if (lo <= k && k < hi) {
						dim = d;
						goto EXIT;
					}
				}

				EXIT:
					bases[dim].push_back(Full[k]);
			}
		}

		return bases;
	};
}


#endif


