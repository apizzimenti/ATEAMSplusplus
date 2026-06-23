
#ifndef ATEAMS_COMMON_H
#define ATEAMS_COMMON_H

#include <SparseRREF/sparse_mat.h>
#include <phat/boundary_matrix.h>
#include <phat/compute_persistence_pairs.h>
#include <phat/representations/default_representations.h>

#include <set>
#include <vector>
#include <map>

#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <string>

namespace ATEAMS {
	/** @cond */
	/** Unsigned long int. */
	typedef unsigned long ulong;

	/** Default type used for indexing in SparseRREF matrices. */
	typedef uint32_t index_t;

	/** Default type used for data storage in SparseRREF matrices. */
	typedef ulong data_t;

	/** SparseRREF matrix. */
	typedef SparseRREF::sparse_mat<data_t, index_t> ZpMatrix;

	/** SparseRREF vector. */
	typedef SparseRREF::sparse_vec<data_t, index_t> ZpVector;

	/** SparseRREF field, \f$\mathrm{GF}(p) \cong \Z/p \Z\f$ for \f$p\f$ prime. */
	typedef SparseRREF::field_t Zp;

	/** SparseRREF index of pivots. */
	typedef SparseRREF::pivot_t<index_t> ZpPivot;

	/** Vector of @ref ATEAMS::ZpMatrix. */
	typedef std::vector<ZpMatrix> ZpMatrices;

	/** Vector of @ref ATEAMS::ZpVector. */
	typedef std::vector<ZpVector> ZpVectors;

	/** Alias of @ref ATEAMS::ZpVectors. */
	typedef std::vector<ZpVector> ZpBasis;

	/** Vector of vector of @ref ATEAMS::ZpBasis. */
	typedef std::vector<ZpBasis> ZpBases;

	/** Flat boundary matrix. */
	typedef std::vector<std::vector<int>> FlatBoundaryMatrix;

	/** PHAT column. */
	typedef std::vector<phat::index> PHATColumn;

	/** PHAT boundary matrix. */
	typedef phat::boundary_matrix<phat::bit_tree_pivot_column> PHATBoundaryMatrix;

	/** PHAT persistence pairs. */
	typedef phat::persistence_pairs PHATPairs;

	/** PHAT persistence algorithm. */
	typedef phat::twist_reduction PHATTwist;
	/** @endcond */
}


#endif
