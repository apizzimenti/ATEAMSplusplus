
#ifndef ATEAMS_TOPOLOGY_BASIS_H
#define ATEAMS_TOPOLOGY_BASIS_H


#include "ATEAMS++/topology/persistence/policies.h"
#include "ATEAMS++/topology/persistence/reduction.h"
#include "SparseRREF/sparse_vec.h"


namespace ATEAMS::topology {
	template <typename RingLike>
	inline SparseBases<RingLike> basis(
		complexes::Complex<RingLike>* complex,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		// Get the full boundary matrix; no need to reindex.
		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseBasis<RingLike> reducedChains(complex->size());
		SparseBases<RingLike> bases(complex->Cells.size());

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

		// Reduction policy.
		auto reductionPolicy = [&resources](
			SparseVector<RingLike>& chain,
			int index,
			int dim
		) {
			return chain.size() > 0 && resources.serial->lookup[helpers::youngestOf<RingLike>(chain)] != 0;
		};

		// Reduction method, since we have to mark which columns we use to
		// eliminate the chain.
		auto reduce = [&R, &Full, &complex, &resources](
			SparseVector<RingLike>& cell,
			SparseVector<RingLike>& reduced
		) {
			SparseVector<RingLike> youngest = Full.rows[resources.serial->lookup[helpers::youngestOf<RingLike>(cell)]];

			typename RingLike::dtype q = *youngest.find(helpers::youngestOf<RingLike>(cell));
			typename RingLike::dtype s = resources.arithmetic->negate[resources.arithmetic->invert[q]];

			// Mark the coefficient used to reduce this chain.
			reduced.push_back(
				resources.serial->lookup[helpers::youngestOf<RingLike>(cell)],
				s
			);

			arithmetic::SparseVectorRescaling<RingLike>(s, youngest, R, resources);
			arithmetic::SparseVectorAddition<RingLike>(cell, youngest, R, resources);
		};

		// Reduce blocks.
		vector<int> endpoints = persistence::policies::traversal::twistFull(complex);

		for (int d=endpoints[0]; d >= endpoints[1]; d--) {
			for (int j=complex->Breaks[d][0]; j < complex->Breaks[d][1]; j++) {
				SparseVector<RingLike>& chain = Full.rows[j];
				SparseVector<RingLike> reduced;

				while (reductionPolicy(chain, j, d)) reduce(chain, reduced);

				if (chain.size() > 0) destructionPolicy(chain, j, d);
				else creationPolicy(j, d);

				reduced.canonicalize();
				reducedChains[j] = reduced;
			}
		}

		// Find essential cycles.
		for (int k : resources.serial->marked) {
			if (resources.serial->lookup[k] == 0) {
				int dim;

				// Determine the dimension of the chain.
				for (int d=0; d < complex->Breaks.size(); d++) {
					int lo = complex->Breaks[d][0];
					int hi = complex->Breaks[d][1];

					if (lo <= k && k < hi) {
						dim = d;
						goto EXIT;
					}
				}

				EXIT:
					// Find the reduced chain, then add a 1 at the index at which
					// the youngest cell was added (i.e. the time the cycle was)
					// created).
					SparseVector<RingLike> r = reducedChains[k];
					r.push_back((INDEX)k, (typename RingLike::dtype)1);
					r.canonicalize();

					bases[dim].push_back(r);
			}
		}

		return bases;
	};


	template <typename RingLike>
	inline SparseBases<RingLike> cobasis(
		complexes::Complex<RingLike>* complex,
		SparseBases<RingLike>& bases,
		Ring* R,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		// Count the total number of basis vectors.
		int tbv = 0;
		for (int i=0; i < bases.size(); i++) tbv += bases[i].size();

		/*
			Create a matrix of the following format:

			_________________________
			|               | 1     |
			| basis vectors |   1   |
			|_______________|_____1_|
			|               |       |
			|  coboundary   |   0   |
			|_______________|_______|

			Then RREF. After solving, the top-left block will be zeroed and
			the bottom-right block will contain cobasis vectors.
		*/
		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseMatrix<RingLike> adjoined(tbv + Full.nrow, Full.ncol + tbv);
		SparseBases<RingLike> cobases(complex->Cells.size(), SparseBasis<RingLike>());

		// Stack the basis vectors on top.
		int tbadded = 0;
		for (int i=0; i < bases.size(); i++) {
			for (int b=0; b < bases[i].size(); b++) {
				SparseVector<RingLike> row = bases[i][b];
				row.push_back((INDEX)Full.ncol+tbadded, (typename RingLike::dtype)1);
				adjoined.rows[tbadded] = row;
				tbadded++;
			}
		}

		// Copy the remaining rows.
		for (int j=0; j < Full.nrow; j++) adjoined[j+tbv] = Full.rows[j];

		printSparseRREFmat<RingLike>(adjoined);
		cout << endl << endl;

		// Row-reduce.
		SparsePivots pivots = arithmetic::SparseMatrixRREF<RingLike>(adjoined, R, resources);
		printSparseRREFmat<RingLike>(adjoined);
		cout << endl << endl;
		arithmetic::SparseMatrixReduce<RingLike>(adjoined, pivots);
		printSparseRREFmat<RingLike>(adjoined);
		cout << endl << endl;

		return cobases;
	}
}


#endif


