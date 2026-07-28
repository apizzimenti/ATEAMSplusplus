
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_POLICIES_H
#define AETAMS_TOPOLOGY_PERSISTENCE_POLICIES_H

#include "ATEAMS++/topology/helpers.h"

namespace ATEAMS::topology::persistence {

	// Reindexing policies
	template <typename RingLike>
	using ReindexingPolicy = function<
		SparseMatrix<RingLike>
		(
			complexes::Complex<RingLike>*,
			vector<int>&,
			arithmetic::ComputeOptions<RingLike>&
		)
	>;

	template <typename RingLike>
	inline SparseMatrix<RingLike> reindexSingle(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		arithmetic::ComputeOptions<RingLike>& options,
		int dimension
	) {
		return helpers::reindexSparseBoundaryMatrix(complex, filtration, dimension, options);
	}

	template <typename RingLike>
	inline SparseMatrix<RingLike> reindexFull(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		return helpers::reindexSparseBoundaryMatrix(complex, filtration, options);
	}

	// Traversal policies
	template <typename RingLike>
	using TraversalPolicy = function<vector<int>(complexes::Complex<RingLike>*)>;

	template <typename RingLike>
	inline vector<int> standardRestrictedTraversalPolicy(
		complexes::Complex<RingLike>* complex,
		int dimension
	) {
		return {
			dimension,
			min(dimension+1, (int)complex->Cells.size()-1)
		};
	}

	template <typename RingLike>
	inline vector<int> standardFullTraversalPolicy(
		complexes::Complex<RingLike>* complex
	) {
		return {
			0,
			(int)complex->Cells.size()-1
		};
	}

	template <typename RingLike>
	inline vector<int> twistRestrictedTraversalPolicy(
		complexes::Complex<RingLike>* complex,
		int dimension
	) {
		return {
			min(dimension+1, (int)complex->Cells.size()-1),
			dimension
		};
	}

	template <typename RingLike>
	inline vector<int> twistFullTraversalPolicy(
		complexes::Complex<RingLike>* complex
	) {
		return {
			(int)complex->Cells.size()-1,
			0
		};
	}


	// Reduction policies
	template <typename RingLike>
	using ReductionPolicy = function<bool(SparseVector<RingLike>&,vector<int>&,int,int)>;


	template <typename RingLike>
	inline bool standardReductionPolicy(
		SparseVector<RingLike>& cell,
		vector<int>& lookup,
		int cellIndex,
		int dim
	) {
		return (cell.size() > 0) && (lookup[helpers::youngestOf<RingLike>(cell)] != 0);
	}

	template <typename RingLike>
	inline bool JITReductionPolicy(
		SparseVector<RingLike>& cell,
		vector<int>& lookup,
		int cellIndex,
		int dim,
		vector<bool>& zeroed
	) {
		return (
			cell.size() > 0 &&
			!zeroed[cellIndex] &&
			lookup[helpers::youngestOf<RingLike>(cell)] != 0
		);
	}





	// Creation policies
	using CreationPolicy = function<void(int,int)>;

	template <typename RingLike>
	inline void standardCreationPolicy(
		int markedIndex,
		int dim,
		set<int>& marked
	) {
		marked.insert(markedIndex);
	}


	template <typename RingLike>
	inline void standardParallelCreationPolicy(
		int markedIndex,
		int dim,
		set<int>& marked,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		options.parallel->marked[dim].insert(markedIndex);
	}





	// Destruction policies
	template <typename RingLike>
	using DestructionPolicy = function<void(SparseVector<RingLike>&,int,int)>;

	template <typename RingLike>
	inline void standardDestructionPolicy(
		SparseVector<RingLike>& cell,
		int markedIndex,
		int dim,
		vector<int>& lookup
	) {
		lookup[helpers::youngestOf<RingLike>(cell)] = markedIndex;
	}

	template <typename RingLike>
	inline void twistDestructionPolicy(
		SparseVector<RingLike>& cell,
		int markedIndex,
		int dim,
		vector<int>& lookup,
		SparseMatrix<RingLike>& Full
	) {
		lookup[helpers::youngestOf<RingLike>(cell)] = markedIndex;
		Full.rows[helpers::youngestOf<RingLike>(cell)].zero();
	};

	template <typename RingLike>
	inline void JITDestructionPolicy(
		SparseVector<RingLike>& cell,
		int markedIndex,
		int dim,
		vector<int>& lookup,
		vector<bool>& zeroed
	) {
		lookup[helpers::youngestOf<RingLike>(cell)] = markedIndex;
		zeroed[markedIndex] = true;
	}



	// Reporting policies
	template <typename RingLike>
	using ReportingPolicy = function<vector<int>(complexes::Complex<RingLike>*,vector<int>&,set<int>&)>;

	template <typename RingLike>
	inline vector<int> standardRestrictedReportingPolicy(
		complexes::Complex<RingLike>* complex,
		vector<int>& lookup,
		set<int>& marked,
		int dimension
	) {
		int low = complex->Breaks[dimension][0];
		int high = complex->Breaks[dimension][1];
		
		vector<int> essential;

		for (int k : marked) {
			if (lookup[k] == 0 && (low <= k && k < high)) essential.push_back(k);
		}

		return essential;
	}


	template <typename RingLike>
	inline vector<int> standardFullReportingPolicy(
		complexes::Complex<RingLike>* complex,
		vector<int>& lookup,
		set<int>& marked
	) {
		
		vector<int> essential;
		for (int k : marked) {
			if (lookup[k] == 0) essential.push_back(k);
		}

		return essential;
	}
}

#endif
