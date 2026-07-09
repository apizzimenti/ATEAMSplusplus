
#ifndef ATEAMS_STATISTICS_OBSERVABLES_H
#define ATEAMS_STATISTICS_OBSERVABLES_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/complexes/Complex.h"

#include <algorithm>
#include <SparseRREF/sparse_mat.h>

namespace ATEAMS {
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
		template <typename T>
		inline int energy(
			ATEAMS::complexes::Complex<T>* complex,
			ATEAMS::SparseVector<T> f,
			ATEAMS::Field F,
			int dimension
		) {
			ATEAMS::SparseVector<T> coefficients = SparseRREF::sparse_mat_dot_sparse_vec<T,INDEX>(
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
		 * 
		 * @returns Number of satisfied cells.
		 */
		template <typename T>
		inline int occupation(
			ATEAMS::complexes::Complex<T>* complex,
			ATEAMS::SparseVector<T> f,
			ATEAMS::Field F,
			int dimension
		) {
			ATEAMS::SparseVector<T> coefficients = SparseRREF::sparse_mat_dot_sparse_vec<T,INDEX>(
				complex->Coboundary.Matrices[dimension], f, F
			);
			coefficients.compress();

			return complex->Cells[dimension]-coefficients.size();
		}

		/**
		 * @brief Computes the indices of satisfied \f$d\f$-cells in the complex
		 * with respect to \f$f_t\f$: that is, the indices of \f$d\f$-cells \f$x\f$
		 * for which \f$(\delta^{d-1}f_t)(x) = 0\f$.
		 * 
		 * @param complex (Pointer to) a Complex.
		 * @param f \f$(d-1)\f$-cochain.
		 * @param F Finite field.
		 * @param dimension Dimension of the cells over which we're sampling.
		 * 
		 * @returns Indices of satisfied cells.
		 */
		template <typename T>
		inline std::vector<int> satisfied(
			ATEAMS::complexes::Complex<T>* complex,
			ATEAMS::SparseVector<T> f,
			ATEAMS::Field F,
			int dimension
		) {
			// Perform the matrix multiplication.
			ATEAMS::SparseVector<T> coefficients = SparseRREF::sparse_mat_dot_sparse_vec<T,INDEX>(
				complex->Coboundary.Matrices[dimension], f, F
			);
			coefficients.compress();

			// Create a vector containing the indices of nonzero entries.
			std::vector<int> unsatisfied(coefficients.size());
			for (int i=0; i < unsatisfied.size(); i++) unsatisfied[i] = (int)coefficients(i);

			// Create a vector containing the indices of *all* entries, then take
			// the difference.
			std::vector<int> totality(complex->Cells[dimension]);
			std::iota(totality.begin(), totality.end(), 0);

			std::vector<int> satisfied;
			std::set_difference(
				totality.begin(), totality.end(),
				unsatisfied.begin(), unsatisfied.end(),
				inserter(satisfied, satisfied.begin())
			);

			return satisfied;
		}

		/**
		 * @brief Computes the indices of unsatisfied \f$d\f$-cells in the complex
		 * with respect to \f$f_t\f$: that is, the indices of \f$d\f$-cells \f$x\f$
		 * for which \f$(\delta^{d-1}f_t)(x) \neq 0\f$.
		 * 
		 * @param complex (Pointer to) a Complex.
		 * @param f \f$(d-1)\f$-cochain.
		 * @param F Finite field.
		 * @param dimension Dimension of the cells over which we're sampling.
		 * 
		 * @returns Indices of unsatisfied cells.
		 */
		template <typename T>
		inline std::vector<int> unsatisfied(
			ATEAMS::complexes::Complex<T>* complex,
			ATEAMS::SparseVector<T> f,
			ATEAMS::Field F,
			int dimension
		) {
			// Perform the matrix multiplication.
			ATEAMS::SparseVector<T> coefficients = SparseRREF::sparse_mat_dot_sparse_vec<T,INDEX>(
				complex->Coboundary.Matrices[dimension], f, F
			);
			coefficients.compress();

			// Create a vector containing the indices of nonzero entries.
			std::vector<int> unsatisfied(coefficients.size());
			for (int i=0; i < unsatisfied.size(); i++) unsatisfied[i] = (int)coefficients(i);

			return unsatisfied;
		}
	}
}

#endif
