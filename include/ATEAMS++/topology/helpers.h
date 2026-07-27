
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_POLICIES_H
#define AETAMS_TOPOLOGY_PERSISTENCE_POLICIES_H

#include <set>
#include <vector>

/** @brief Helper functionality. */
namespace ATEAMS::topology::helpers {

	template <typename RingLike>
	inline INDEX youngestOf(SparseVector<RingLike> cell) {
		return cell(cell.size()-1);
	}

	template <typename RingLike>
	inline SparseMatrix<RingLike> reindexSparseBoundaryMatrix(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		// Construct an index mapping.
		vector<int> remapping(filtration.size());
		
		#pragma omp parallel for shared(filtration, remapping) if(options.parallel->enabled)
		for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

		// Duplicate the matrix.
		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseMatrix<RingLike> Reindexed(Full.nrow, Full.ncol);

		#pragma omp parallel for shared(filtration, remapping)
		for (int t=0; t < Full.nrow; t++) {
			SparseVector<RingLike>& row = Reindexed.rows[filtration[t]];
			SparseVector<RingLike> orow = Full.rows[t];

			for (int i=0; i < orow.size(); i++) {
				row.push_back(
					(INDEX)remapping[orow(i)],
					(typename RingLike::dtype)orow[i]
				);
			}
		}

		Reindexed.compress();
		return Reindexed;
	}

	
	template <typename RingLike>
	inline SparseMatrix<RingLike> reindexSparseBoundaryMatrix(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		int dimension,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		// Construct an index mapping.
		vector<int> remapping(filtration.size());
		
		#pragma omp parallel for shared(filtration, remapping) if(options.parallel->enabled)
		for (int t=0; t < filtration.size(); t++) remapping[filtration[t]] = t;

		SparseMatrix<RingLike> Full = complex->Coboundary.Full;
		SparseMatrix<RingLike> Reindexed(Full.nrow, Full.ncol);

		int startDimension = complex->Breaks[dimension][0];
		int stopDimension = complex->Breaks[dimension][1];

		#pragma omp parallel for shared(Full, Reindexed) if (options.parallel->enabled)
		for (int t=0; t < Full.nrow; t++) {
			if ((startDimension <= t) && (t < stopDimension)) {
				Reindexed.rows[t] = Full.rows[filtration[t]];
			} else {
				SparseVector<RingLike>& row = Reindexed.rows[t];
				SparseVector<RingLike> orow = Full.rows[t];

				for (int i=0; i < orow.size(); i++) {
					row.push_back(
						(INDEX)remapping[orow(i)],
						(typename RingLike::dtype)orow[i]
					);
				}
			}
		}

		Reindexed.compress();
		return Reindexed;
	}
}




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
	inline SparseMatrix<RingLike> reindexWhole(
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
			min(dimension+1, (int)complex->Cells.size())
		};
	}

	template <typename RingLike>
	inline vector<int> standardFullTraversalPolicy(
		complexes::Complex<RingLike>* complex
	) {
		return { 0, (int)complex->Cells.size() };
	}

	template <typename RingLike>
	inline vector<int> twistRestrictedTraversalPolicy(
		complexes::Complex<RingLike>* complex,
		int dimension
	) {
		return {
			min(dimension+1, (int)complex->Cells.size()),
			dimension
		};
	}

	template <typename RingLike>
	inline vector<int> twistFullTraversalPolicy(
		complexes::Complex<RingLike>* complex
	) {
		return { (int)complex->Cells.size(), 0 };
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
		vector<bool> zeroed
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