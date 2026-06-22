
#ifndef ATEAMS_STATISTICS_OBSERVABLES_H
#define ATEAMS_STATISTICS_OBSERVABLES_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/complexes/complexes.h"

#include "libraries/SparseRREF/sparse_mat.h"

namespace ATEAMS {
	/** @brief Statistical support. */
	namespace statistics {

		/**
		 * @brief Computes the total energy of the given configuration.
		 * 
		 * @param complex (Pointer to) a Complex.
		 * @param f \f$(d-1)\f$-cochain.
		 * @param F Finite field.
		 * @param dimension Dimension of the cells over which we're sampling.
		 * 
		 * @returns An integer giving the (negative) number of \f$(d-1)\f$-cells
		 * over which \f$f_t\f$ does not vanish.
		 */
		inline int energy(
			ATEAMS::complexes::Complex* complex,
			ATEAMS::ZpVector f,
			ATEAMS::Zp F,
			int dimension
		) {
			ATEAMS::ZpVector coefficients = SparseRREF::sparse_mat_dot_sparse_vec<data_t,index_t>(
				complex->Coboundary.Matrices[dimension], f, F
			);
			coefficients.compress();

			return -coefficients.size();
		}

		/**
		 * @brief Computes the number of satisfied \f$d\f$-cells in the complex
		 * with respect to \f$f_t\f$: that is, the number of \f$d\f$-cells \f$x\f$
		 * for which \f$(\delta^{d-1}f_t)(x) = 0\f$.
		 * 
		 * @param complex (Pointer to) a Complex.
		 * @param f \f$(d-1)\f$-cochain.
		 * @param F Finite field.
		 * @param dimension Dimension of the cells over which we're sampling.
		 */
		inline int occupation(
			ATEAMS::complexes::Complex* complex,
			ATEAMS::ZpVector f,
			ATEAMS::Zp F,
			int dimension
		) {
			ATEAMS::ZpVector coefficients = SparseRREF::sparse_mat_dot_sparse_vec<data_t,index_t>(
				complex->Coboundary.Matrices[dimension], f, F
			);
			coefficients.compress();

			return complex->Cells[dimension]-coefficients.size();
		}
	};
}

#endif
