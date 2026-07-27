
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_T
#define ATEAMS_TOPOLOGY_PERSISTENCE_T

#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_H
#error __FILE__ should only be included from topology/persistence.h.
#endif

#include "ATEAMS++/topology/persistence.h"
#include "ATEAMS++/arithmetic/arithmetic.h"

#include <phat/compute_persistence_pairs.h>
#include <phat/helpers/misc.h>

using namespace SparseRREF;
using namespace std;

namespace ATEAMS::topology::persistence {

	/** @cond */
	template <typename RingLike>
	inline void reduceChain(
		SparseMatrix<RingLike>& Full,
		vector<int>& lookup,
		SparseVector<RingLike>& cell,
		Ring* R
	) {
		// Get the youngest chain (column) sharing a face with `cell` and zero
		// out the pivot row.
		SparseVector<RingLike> youngest = Full.rows[lookup[helpers::youngestOf<RingLike>(cell)]];

		typename RingLike::dtype q = *youngest.find(helpers::youngestOf<RingLike>(cell));
		typename RingLike::dtype s = scalar_neg(scalar_inv(q, R->ring), R->ring);

		arithmetic::SparseVectorRescaling<RingLike>(s, youngest, R);
		arithmetic::SparseVectorAddition<RingLike>(cell, youngest, R);
	}


	template <typename RingLike>
	inline void reduceBlock(
		SparseMatrix<RingLike>& Full,
		vector<int>& endpoints,
		vector<int>& lookup,
		int dim,
		Ring *R,
		ReductionPolicy<RingLike> reductionPolicy,
		CreationPolicy creationPolicy,
		DestructionPolicy<RingLike> destructionPolicy,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		for (int j=endpoints[0]; j < endpoints[1]; j++) {
			SparseVector<RingLike>& cell = Full.rows[j];

			while (reductionPolicy(cell, lookup, j, dim)) reduceChain<RingLike>(Full, lookup, cell, R);

			if (cell.size() > 0) destructionPolicy(cell, j, dim);
			else creationPolicy(j, dim);
		}
	}

	/** @endcond */
}


namespace ATEAMS::topology::persistence {

	template <typename RingLike>
	inline vector<int> standard(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options,
		ReindexingPolicy<RingLike> reindexingPolicy,
		TraversalPolicy<RingLike> traversalPolicy,
		ReportingPolicy<RingLike> reportingPolicy
	) {
		// Determine the endpoints and reindex the boundary matrix accordingly.
		SparseMatrix<RingLike> Full = reindexingPolicy(complex, filtration, options);

		// Track which column is to be added next; track which ones are marked.
		vector<int> youngestChainLookup(complex->size(), 0);
		set<int> marked;

		// Cycle creation policy.
		auto creationPolicy = std::bind(
			standardCreationPolicy<RingLike>,	// using the standard creation policy
			placeholders::_1,					// placeholder for `markedIndex`
			placeholders::_2,					// placeholder for `dim`
			std::ref(marked)					// a reference to `marked`.
		);

		// Cycle destruction policy.
		auto destructionPolicy = std::bind(
			standardDestructionPolicy<RingLike>,	// using the standard destruction policy
			placeholders::_1,						// placeholder for `cell`
			placeholders::_2,						// placeholder for `markedIndex`
			placeholders::_3,						// placeholder for `dim`,
			std::ref(youngestChainLookup)			// reference to `youngestChainLookup`.
		);

		// Reduce the blocks of the matrix.
		vector<int> endpoints = traversalPolicy(complex);

		for (int d=endpoints[0]; d <= endpoints[1]; d++) {
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				youngestChainLookup,
				d,
				R,
				standardReductionPolicy<RingLike>,
				creationPolicy,
				destructionPolicy,
				options
			);
		}

		// Find essential cycles.
		return reportingPolicy(
			complex,
			youngestChainLookup,
			marked
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
			standardRestrictedTraversalPolicy<RingLike>,	// standard restricted traversal policy, since we're in a specific range
			placeholders::_1,								// placeholder for `complex`
			dimension										// autofill the `dimension` parameter.
		);

		auto reportingPolicy = std::bind(
			standardRestrictedReportingPolicy<RingLike>,	// again restricted, since we're in a range
			placeholders::_1,
			placeholders::_2,
			placeholders::_3,
			dimension
		);

		auto reindexingPolicy = std::bind(
			reindexSingle<RingLike>,
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
	inline vector<int> twist(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options,
		ReindexingPolicy<RingLike> reindexingPolicy,
		TraversalPolicy<RingLike> traversalPolicy,
		ReportingPolicy<RingLike> reportingPolicy
	) {
		// Determine the endpoints and reindex the boundary matrix accordingly.
		SparseMatrix<RingLike> Full = reindexingPolicy(complex, filtration, options);

		// Track which column is to be added next; track which ones are marked.
		vector<int> youngestChainLookup(complex->size(), 0);
		set<int> marked;

		// Cycle creation policy.
		auto creationPolicy = std::bind(
			standardCreationPolicy<RingLike>,	// using the standard creation policy
			placeholders::_1,					// placeholder for `markedIndex`
			placeholders::_2,					// placeholder for `dim`
			std::ref(marked)					// a reference to `marked`.
		);

		// Cycle destruction policy.
		auto destructionPolicy = std::bind(
			twistDestructionPolicy<RingLike>,	// using the twist destruction policy
			placeholders::_1,					// placeholder for `cell`
			placeholders::_2,					// placeholer for `markedIndex`
			placeholders::_3,					// placeholder for `dim`
			std::ref(youngestChainLookup),		// reference to `youngestChainLookup`
			std::ref(Full)						// reference to `Full`, for clearing
		);

		// Reduce the blocks of the matrix.
		vector<int> endpoints = traversalPolicy(complex);

		for (int d=endpoints[0]; d >= endpoints[1]; d--) {
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				youngestChainLookup,
				d,
				R,
				standardReductionPolicy<RingLike>,
				creationPolicy,
				destructionPolicy,
				options
			);
		}

		// Find essential cycles.
		return reportingPolicy(
			complex,
			youngestChainLookup,
			marked
		);
	};

	template <typename RingLike>
	inline vector<int> twist(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		auto traversalPolicy = std::bind(
			twistRestrictedTraversalPolicy<RingLike>,		// standard restricted traversal policy, since we're in a specific range
			placeholders::_1,								// placeholder for `complex`
			dimension										// autofill the `dimension` parameter.
		);

		auto reportingPolicy = std::bind(
			standardRestrictedReportingPolicy<RingLike>,	// again restricted, since we're in a range
			placeholders::_1,
			placeholders::_2,
			placeholders::_3,
			dimension
		);

		auto reindexingPolicy = std::bind(
			reindexSingle<RingLike>,
			placeholders::_1,
			placeholders::_2,
			placeholders::_3,
			dimension
		);

		return twist<RingLike>(
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
	inline vector<int> standardParallel(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options,
		ReindexingPolicy<RingLike> reindexingPolicy,
		TraversalPolicy<RingLike> traversalPolicy,
		ReportingPolicy<RingLike> reportingPolicy
	) {
		// Determine the endpoints and reindex the boundary matrix accordingly.
		SparseMatrix<RingLike> Full = reindexingPolicy(complex, filtration, options);

		// Track which column is to be added next; track which ones are marked.
		vector<int> youngestChainLookup(complex->size(), 0);
		set<int> marked;

		// Flush data structures.
		options.parallel->flush();

		// Cycle creation policy.
		auto creationPolicy = std::bind(
			standardParallelCreationPolicy<RingLike>,	// using the standard parallel creation policy
			placeholders::_1,							// placeholder for `markedIndex`
			placeholders::_2,							// placeholder for `dim`
			std::ref(marked),							// a dummy reference to `marked`
			std::ref(options)							// a reference to `options`, for marking in parallel.
		);

		// Cycle destruction policy.
		auto destructionPolicy = std::bind(
			standardDestructionPolicy<RingLike>,	// using the standard destruction policy
			placeholders::_1,						// placeholder for `cell`
			placeholders::_2,						// placeholder for `markedIndex`
			placeholders::_3,						// placeholder for `dim`,
			std::ref(youngestChainLookup)			// reference to `youngestChainLookup`.
		);

		// Reduce the blocks in parallel, since they are independent of one another.
		vector<int> endpoints = traversalPolicy(complex);

		#pragma omp parallel for firstprivate(Full) shared(youngestChainLookup)
		for (int d=endpoints[0]; d <= endpoints[1]; d++) {
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				youngestChainLookup,
				d,
				R,
				standardReductionPolicy<RingLike>,
				creationPolicy,
				destructionPolicy,
				options
			);
		}

		// Re-constitute the marked columns.
		for (int d=endpoints[0]; d <= endpoints[1]; d++) {
			for (auto& k : options.parallel->marked[d]) marked.insert(k);
		}

		// Find essential cycles.
		return reportingPolicy(
			complex,
			youngestChainLookup,
			marked
		);
	}


	template <typename RingLike>
	inline vector<int> standardParallel(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		auto traversalPolicy = std::bind(
			standardRestrictedTraversalPolicy<RingLike>,		// standard restricted traversal policy, since we're in a specific range
			placeholders::_1,								// placeholder for `complex`
			dimension										// autofill the `dimension` parameter.
		);

		auto reportingPolicy = std::bind(
			standardRestrictedReportingPolicy<RingLike>,	// again restricted, since we're in a range
			placeholders::_1,
			placeholders::_2,
			placeholders::_3,
			dimension
		);

		auto reindexingPolicy = std::bind(
			reindexSingle<RingLike>,
			placeholders::_1,
			placeholders::_2,
			placeholders::_3,
			dimension
		);

		return standardParallel<RingLike>(
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
		arithmetic::ComputeOptions<RingLike>& options,
		ReindexingPolicy<RingLike> reindexingPolicy,
		TraversalPolicy<RingLike> traversalPolicy,
		ReportingPolicy<RingLike> reportingPolicy
	) {
		// Determine the endpoints and reindex the boundary matrix accordingly.
		SparseMatrix<RingLike> Full = reindexingPolicy(complex, filtration, options);

		// Track which column is to be added next; track which ones are zeroed;
		// track which ones are marked.
		vector<int> youngestChainLookup(complex->size(), 0);
		vector<bool> zeroed(complex->size(), 0);
		set<int> marked;

		// Cycle creation policy.
		auto creationPolicy = std::bind(
			standardParallelCreationPolicy<RingLike>,	// using the standard parallel creation policy
			placeholders::_1,							// placeholder for `markedIndex`
			placeholders::_2,							// placeholder for `dim`
			std::ref(marked),							// a dummy reference to `marked`
			std::ref(options)							// a reference to `options`, for marking in parallel.
		);

		// Cycle destruction policy.
		auto destructionPolicy = std::bind(
			JITDestructionPolicy<RingLike>,	// using the JIT destruction policy
			placeholders::_1,				// placeholder for `cell`
			placeholders::_2,				// placeholder for `markedIndex`
			placeholders::_3,				// placeholder for `dim`,
			std::ref(youngestChainLookup),	// reference to `youngestChainLookup`,
			std::ref(zeroed)				// reference to `zeroed`.
		);

		// Reduction policy.
		auto reductionPolicy = std::bind(
			JITReductionPolicy<RingLike>,	// using the JIT reduction policy
			placeholders::_1,				// placeholder for `cell`						
			placeholders::_2,				// placeholder for `lookup`
			placeholders::_3,				// placeholder for `cellIndex`
			placeholders::_4,				// placeholder for `dim`
			std::ref(zeroed)				// reference to `zeroed`.
		);

		// Reduce the blocks in parallel, since they are independent of one another.
		vector<int> endpoints = traversalPolicy(complex);

		#pragma omp parallel for firstprivate(Full) shared(youngestChainLookup,zeroed)
		for (int d=endpoints[0]; d <= endpoints[1]; d++) {
			reduceBlock<RingLike>(
				Full,
				complex->Breaks[d],
				youngestChainLookup,
				d,
				R,
				reductionPolicy,
				creationPolicy,
				destructionPolicy,
				options
			);
		}

		// Re-constitute the marked columns.
		for (int d=endpoints[0]; d <= endpoints[1]; d++) {
			for (auto& k : options.parallel->marked[d]) marked.insert(k);
		}

		// Find essential cycles.
		return reportingPolicy(
			complex,
			youngestChainLookup,
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
			standardRestrictedTraversalPolicy<RingLike>,		// standard restricted traversal policy, since we're in a specific range
			placeholders::_1,								// placeholder for `complex`
			dimension										// autofill the `dimension` parameter.
		);

		auto reportingPolicy = std::bind(
			standardRestrictedReportingPolicy<RingLike>,	// again restricted, since we're in a range
			placeholders::_1,
			placeholders::_2,
			placeholders::_3,
			dimension
		);

		auto reindexingPolicy = std::bind(
			reindexSingle<RingLike>,
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
	inline vector<int> PHAT(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration
	) {
		// The filtration specifies the order in which we add the cells of all
		// dimensions. Create a map that specifies to which position each cell was
		// moved. For example, if the filtration has
		//
		//				[... 12, 9, 10, 19, ...]
		//
		// that was originally
		//				[... 9, 10, 11, 12, ...]
		//
		// then the mapping should have entries
		//
		//				{... 12: 9, 9: 10, 10: 11, 19: 12, ...}

		// TODO optimize, this takes way too long
		vector<int> remapping(filtration.size(), 0);
		for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

		PHATBoundaryMatrix boundary;
		PHATColumn column;

		vector<int> ground;
		int dim;

		boundary.set_num_cols(complex->Boundary.Flat.size());

		// this also takes an insane amount of time; no need to go all the way
		// up.
		for (int t=0; t < complex->Boundary.Flat.size(); t++) {
			ground = complex->Boundary.Flat[filtration[t]];

			for (int j=0; j < ground.size(); j++) {
				column.push_back(remapping[ground[j]]);
			}
			// Check whether we're in the dimension of the cells on which we're
			// percolating; if so, add the cells in the order specified in the
			// filtration. If we're in the percolation dimension + 1, re-index the
			// (d+1)-dim cells to reflect the ordering of the d-dimensional cells.
			// Otherwise, do nothing.
			// if ((startDimension <= t) && (t < stopDimension)) {
			// 	for (int j=0; j < faces; j++) column.push_back(complex->Boundary.Flat[filtration[t]][j]);
			// } else {
			// 	for (int j=0; j < faces; j++) column.push_back(remapping[complex->Boundary.Flat[t][j]]);
			// }

			// for (int j=0; j < faces; j++) {
			// 	column.push_back()
			// }

			std::sort(column.begin(), column.end());
		
			// Determine the dimension of the cell.
			for (int d=0; d < complex->Breaks.size(); d++) {
				int lo = complex->Breaks[d][0];
				int hi = complex->Breaks[d][1];

				if (lo <= filtration[t] && filtration[t] < hi) {
					dim = d;
					break;
				}
			}

			boundary.set_dim(t, dim);
			boundary.set_col(t, column);
			column.clear();
		}

		// Compute the persistence pairs and populate a Vector to return to the
		// user. We know when the giant(/essential) cycles are born because they
		// do not appear in the list of births or deaths.
		PHATPairs pairs;
		set<int> all;
		phat::compute_persistence_pairs<PHATTwist>(pairs, boundary);
		pairs.sort();

		for (phat::index i=0; i < pairs.get_num_pairs(); i++) {
			all.insert(pairs.get_pair(i).first);
			all.insert(pairs.get_pair(i).second);
		}

		vector<int> essential;

		for (int t=0; t < filtration.size(); t++) {
			if (!all.contains(filtration[t])) essential.push_back(filtration[t]);
		}

		return essential;
	}



	template <typename RingLike>
	inline vector<int> PHAT(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		int dimension
	) {
		// No sense in doing a ton of policy work for this one.
		vector<int> essential = PHAT<RingLike>(complex, filtration);
		int low = complex->Breaks[dimension][0], high = complex->Breaks[dimension][1];
		std::erase_if(essential, [filtration,low,high](int t) { return !(low <= filtration[t] && filtration[t] < high); });

		return essential;
	}


	template <typename RingLike>
	std::vector<int> persistence(
		complexes::Complex<RingLike>* complex,
		std::vector<int>& filtration,
		Ring* R,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		vector<int> essential;

		if (R->characteristic == 2) essential = PHAT<RingLike>(complex, filtration, dimension);
		else {
			if (options.parallel->enabled) essential = JIT<RingLike>(complex, filtration, R, dimension, options);
			else essential = twist<RingLike>(complex, filtration, R, dimension, options);
		}

		std::sort(essential.begin(), essential.end());
		return essential;
	};
}

#endif
