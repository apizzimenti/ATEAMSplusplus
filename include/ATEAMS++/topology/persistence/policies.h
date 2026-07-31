
#ifndef ATEAMS_TOPOLOGY_PERSISTENCE_POLICIES_H
#define AETAMS_TOPOLOGY_PERSISTENCE_POLICIES_H

#include "ATEAMS++/topology/helpers.h"

namespace ATEAMS::topology::persistence::policies::reindexing {

	/**
	 * @brief Reindexes only the columns corresponding to the in-focus dimension.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix. For example, if we are doing
	 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
	 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
	 * @param resources Computing resources.
	 * @param dimension The percolation dimension.
	 * 
	 * @returns A reindexed sparse boundary matrix.
	 */
	template <typename RingLike>
	inline SparseMatrix<RingLike> single(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		arithmetic::ComputeResources<RingLike>& resources,
		int dimension
	) {
		return helpers::reindexSparseBoundaryMatrix(complex, filtration, dimension, resources);
	}

	/**
	 * @brief Reindexes an entire sparse boundary matrix.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
	 * @param filtration A vector that specifies the order in which to add
	 * 	the cells in the flat boundary matrix. For example, if we are doing
	 * 	1-dimensional percolation in 2-d, then we switch up the ordering of
	 * 	the 1-d cells in the flat boundary matrix, but leave everything else.
	 * @param resources Computing resources.
	 * 
	 * @returns A reindexed sparse boundary matrix.
	 */
	template <typename RingLike>
	inline SparseMatrix<RingLike> full(
		complexes::Complex<RingLike>* complex,
		vector<int>& filtration,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		return helpers::reindexSparseBoundaryMatrix(complex, filtration, resources);
	}
}



namespace ATEAMS::topology::persistence::policies::traversal {
	/**
	 * @brief For a given dimension \f$d\f$, gives initial and terminal values
	 * \f$[d, \min(d+1, \max(X)-1)]\f$, where \f$\max(X)\f$ is the maximal dimension
	 * of the complex \f$X\f$.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
	 * @param dimension In-focus dimension.
	 * 
	 * @returns Initial and terminal values \f$[d, \min(d+1, \max(X)-1)]\f$,
	 * where \f$\max(X)\f$ is the maximal dimension of the complex \f$X\f$.
	 */
	template <typename RingLike>
	inline vector<int> standardRestricted(
		complexes::Complex<RingLike>* complex,
		int dimension
	) {
		return {
			dimension,
			min(dimension+1, (int)complex->Cells.size()-1)
		};
	}


	/**
	 * @brief For a given dimension \f$d\f$, gives initial and terminal values
	 * \f$[0, \max(X)-1]\f$, where \f$\max(X)\f$ is the maximal dimension
	 * of the complex \f$X\f$.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
	 * 
	 * @returns Initial and terminal values \f$[0, \max(X)-1]\f$.
	 */
	template <typename RingLike>
	inline vector<int> standardFull(
		complexes::Complex<RingLike>* complex
	) {
		return {
			0,
			(int)complex->Cells.size()-1
		};
	}

	/**
	 * @brief For a given dimension \f$d\f$, gives initial and terminal values
	 * \f$[\min(d+1, \max(X)-1), d]\f$, where \f$\max(X)\f$ is the maximal dimension
	 * of the complex \f$X\f$ — that is, traverses dimension in reverse order.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
	 * @param dimension In-focus dimension.
	 * 
	 * @returns Initial and terminal values \f$[\min(d+1, \max(X)-1), d]\f$,
	 * where \f$\max(X)\f$ is the maximal dimension of the complex \f$X\f$.
	 */
	template <typename RingLike>
	inline vector<int> twistRestricted(
		complexes::Complex<RingLike>* complex,
		int dimension
	) {
		return {
			min(dimension+1, (int)complex->Cells.size()-1),
			dimension
		};
	}

	/**
	 * @brief For a given dimension \f$d\f$, gives initial and terminal values
	 * \f$[\max(X)-1, 0]\f$, where \f$\max(X)\f$ is the maximal dimension
	 * of the complex \f$X\f$ — that is, traverses dimension in reverse order.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
	 * 
	 * @returns Initial and terminal values \f$[\max(X)-1, 0]\f$,
	 * where \f$\max(X)\f$ is the maximal dimension of the complex \f$X\f$.
	 */
	template <typename RingLike>
	inline vector<int> twistFull(
		complexes::Complex<RingLike>* complex
	) {
		return {
			(int)complex->Cells.size()-1,
			0
		};
	}
}


namespace ATEAMS::topology::persistence::policies::reduction {
	/**
	 * @brief Standard reduction policy: while @p chain has nonzero boundary
	 * and shares its youngest face with another chain, this chain can still be
	 * reduced.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param chain Chain to be checked.
	 * @param lookup Lookup table to find the column corresponding to @p chain's
	 * youngest face.
	 * @param index Column index of the chain.
	 * @param dim Dimension.
	 * 
	 * @returns Whether the chain (column) can be reduced further.
	 */
	template <typename RingLike>
	inline bool standard(
		SparseVector<RingLike>& chain,
		vector<int>& lookup,
		int index,
		int dim
	) {
		return (chain.size() > 0) && (lookup[helpers::youngestOf<RingLike>(chain)] != 0);
	}

	/**
	 * @brief A parallel-friendly implementation of the clearing optimization:
	 * while @p chain has nonzero boundary and shares its youngest face with
	 * another chain, this chain can still be reduced. @p zeroed tracks which
	 * chains have been zeroed by other processes.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param chain Chain to be checked.
	 * @param lookup Lookup table to find the column corresponding to @p chain's
	 * youngest face.
	 * @param index Column index of the chain.
	 * @param dim Dimension.
	 * @param zeroed Boolean vector storing whether a different concurrent
	 * process has cleared this chain.
	 * 
	 * @returns Whether the chain (column) can be reduced further.
	 */
	template <typename RingLike>
	inline bool JIT(
		SparseVector<RingLike>& chain,
		vector<int>& lookup,
		int index,
		int dim,
		vector<bool>& zeroed
	) {
		return (
			chain.size() > 0 &&
			!zeroed[index] &&
			lookup[helpers::youngestOf<RingLike>(chain)] != 0
		);
	}
}


namespace ATEAMS::topology::persistence::policies::creation {
	/**
	 * @brief Standard cycle-creation policy (marking).
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param index Column index of the chain inducing a cycle.
	 * @param dim Dimension.
	 * @param marked Set keeping track of marked chains.
	 */
	template <typename RingLike>
	inline void standard(
		int index,
		int dim,
		set<int>& marked
	) {
		marked.insert(index);
	}

	/**
	 * @brief Standard cycle-creation policy (marking).
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param index Column index of the chain inducing a cycle.
	 * @param dim Dimension.
	 * @param resources Computing resources.
	 */
	template <typename RingLike>
	inline void parallel(
		int index,
		int dim,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		resources.parallel->marked[dim].insert(index);
	}
}



namespace ATEAMS::topology::persistence::policies::destruction {
	/**
	 * @brief A policy that specifies how data for marking boundaries (i.e. the
	 * destruction of cycles) is disseminated.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * Policies of this kind take three arguments: a chain (represented as a
	 * @ref ATEAMS::SparseVector), an integer indicating the column index of the
	 * chain, and an integer corresponding to the dimension of the chain. All
	 * other destruction policies forward to this one.
	 */
	template <typename RingLike>
	using DestructionPolicy = function<void(SparseVector<RingLike>&,int,int)>;

	/**
	 * @brief Standard destruction policy: indicates that this chain destroys
	 * the cycle created by its youngest face.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param chain Chain with boundary.
	 * @param index Column index of the chain.
	 * @param dim Dimension of the chain.
	 * @param lookup Lookup table.
	 */
	template <typename RingLike>
	inline void standard(
		SparseVector<RingLike>& chain,
		int index,
		int dim,
		vector<int>& lookup
	) {
		lookup[helpers::youngestOf<RingLike>(chain)] = index;
	}

	/**
	 * @brief Clearing destruction policy: since the youngest face of this chain
	 * induces a cycle (i.e. the boundary of this chain), we can pre-emptively
	 * clear the column corresponding to the youngest face. This prevents the
	 * boundary column from being processed and immediately marks it as a
	 * persistence pair.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param chain Chain with boundary.
	 * @param index Column index of the chain.
	 * @param dim Dimension of the chain.
	 * @param lookup Lookup table.
	 * @param Full Full (co)boundary matrix.
	 */
	template <typename RingLike>
	inline void twist(
		SparseVector<RingLike>& chain,
		int index,
		int dim,
		vector<int>& lookup,
		SparseMatrix<RingLike>& Full
	) {
		lookup[helpers::youngestOf<RingLike>(chain)] = index;
		Full.rows[helpers::youngestOf<RingLike>(chain)].zero();
	};

	
	/**
	 * @brief A parallel-friendly clearing destruction policy: since the youngest
	 * face of this chain induces a cycle (i.e. the boundary of this chain), we
	 * can pre-emptively mark the column corresponding to its youngest face as
	 * cleared. This prevents the column from being unnecessarily processed.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param chain Chain with boundary.
	 * @param index Column index of the chain.
	 * @param dim Dimension of the chain.
	 * @param lookup Lookup table.
	 * @param zeroed Indicates which columns have been zeroed.
	 */
	template <typename RingLike>
	inline void JIT(
		SparseVector<RingLike>& chain,
		int index,
		int dim,
		vector<int>& lookup,
		vector<bool>& zeroed
	) {
		lookup[helpers::youngestOf<RingLike>(chain)] = index;
		zeroed[helpers::youngestOf<RingLike>(chain)] = true;
	}

	/**
	 * @brief A parallel-friendly clearing destruction policy: since the youngest
	 * face of this chain induces a cycle (i.e. the boundary of this chain), we
	 * can pre-emptively mark the column corresponding to its youngest face as
	 * cleared and, since the concurrent threads' memory locations do not overlap,
	 * actually zero the column. This prevents the column from being unnecessarily
	 * processed.
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param chain Chain with boundary.
	 * @param index Column index of the chain.
	 * @param dim Dimension of the chain.
	 * @param lookup Lookup table.
	 * @param zeroed Indicates which columns have been zeroed.
	 * @param Full Full (co)boundary matrix.
	 */
	template <typename RingLike>
	inline void split(
		SparseVector<RingLike>& chain,
		int index,
		int dim,
		vector<int>& lookup,
		vector<bool>& zeroed,
		SparseMatrix<RingLike>& Full
	) {
		lookup[helpers::youngestOf<RingLike>(chain)] = index;
		zeroed[helpers::youngestOf<RingLike>(chain)] = true;
		Full.rows[helpers::youngestOf<RingLike>(chain)].zero();
	}
}






namespace ATEAMS::topology::persistence::policies::reporting {
	/**
	 * @brief Standard reporting policy, restricted to a range: essential(/giant)
	 * cycles are cycles that are not destroyed (i.e. cycles that are not boundaries).
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
	 * @param lookup Lookup table.
	 * @param marked Set of column indices corresponding to cycles.
	 * @param dimension In-focus dimension.
	 * 
	 * @returns Vector of times at which @p dimension-dimensional essential cycles
	 * were created.
	 */
	template <typename RingLike>
	inline vector<int> standardRestricted(
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

	/**
	 * @brief Standard reporting policy, restricted to a range: essential(/giant)
	 * cycles are cycles that are not destroyed (i.e. cycles that are not boundaries).
	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
	 * 
	 * @param complex A (pointer to) a @ref ATEAMS::complexes::Complex.
	 * @param lookup Lookup table.
	 * @param marked Set of column indices corresponding to cycles.
	 * 
	 * @returns Vector of times at which essential cycles were created.
	 */
	template <typename RingLike>
	inline vector<int> standardFull(
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
