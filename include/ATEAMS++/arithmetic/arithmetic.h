
#ifndef ATEAMS_ARITHMETIC_H
#define ATEAMS_ARITHMETIC_H

#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_vec.h>
#include <SparseRREF/scalar.h>

#include "ATEAMS++/arithmetic/options.h"

namespace ATEAMS {
	/**
	 * @brief Sparse, finite-field linear algebra and support routines.
	 */
	namespace arithmetic {

		/**
		 * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over a
		 * @ref ATEAMS::Ring
	 	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) the coefficient ring @ref Q.
		 */
		template <typename RingLike>
		inline void SparseVectorAddition(
			SparseVector<RingLike>& u,
			SparseVector<RingLike>& v,
			Ring* R
		);

		/**
		 * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Z2.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) the coefficient ring @ref Q.
		 */
		template <>
		inline void SparseVectorAddition<Z2>(
			SparseVector<Z2>& u,
			SparseVector<Z2>& v,
			Ring* R
		) {
			sparse_vec_add<INDEX>(u, v, R->ring);
		};

		/**
		 * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Zp.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) the coefficient ring @ref Q.
		 */
		template <>
		inline void SparseVectorAddition<Zp>(
			SparseVector<Zp>& u,
			SparseVector<Zp>& v,
			Ring* R
		) {
			sparse_vec_add<INDEX>(u, v, R->ring);
		};

		/**
		 * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Q.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) the coefficient ring @ref Q.
		 */
		template <>
		inline void SparseVectorAddition<Q>(
			SparseVector<Q>& u,
			SparseVector<Q>& v,
			Ring* R
		) {
			sfmpq_vec_addsub_mul<INDEX,false>(u, v, (Q::dtype)1);
		}

		/**
		 * @brief Sparse (right) matrix multiplication \f$ A \vec x \f$.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param A Matrix.
		 * @param x Vector.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @returns \f$ \vec y = A \vec x \f$.
		 */
		template <typename RingLike>
		inline SparseVector<RingLike> SparseRightMultiplication(
			SparseMatrix<RingLike>& A,
			SparseVector<RingLike>& x,
			Ring* R
		) {
			SparseVector<RingLike> y = sparse_mat_dot_sparse_vec<typename RingLike::dtype,INDEX>(
				A, x, R->ring
			);
			
			return y;
		};


		/**
		 * @brief Sparse vector rescaling \f$ a \vec x \f$.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param a Scalar.
		 * @param x Vector.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 */
		template <typename RingLike>
		inline void SparseVectorRescaling(
			typename RingLike::dtype a,
			SparseVector<RingLike>& x,
			Ring* R
		) {
			sparse_vec_rescale<INDEX,typename RingLike::dtype>(x, a, R->ring);
		};

		/**
		 * @brief Computes the reduced row echelon form \f$ \RREF(A) \f$.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param A Matrix.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * @param options Multithreaded computing options.
		 * 
		 * @returns \f$ \RREF(A) \f$.
		 */
		template <typename RingLike>
		inline SparsePivots SparseMatrixRREF(
			SparseMatrix<RingLike>& A,
			Ring* R,
			ComputeResources<RingLike>& options
		) {
			return SparseRREF::sparse_mat_rref<typename RingLike::dtype,INDEX>(A, R->ring, options.opt);
		};

		/**
		 * @brief Computes \f$ \ker(A) \f$, where \f$A\f$ is in reduced row
		 * 	echelon form.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param A Matrix.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * @param pivots Pivot indices of \f$ A \f$.
		 * @param options Multithreaded computing options.
		 * 
		 * @returns \f$ \ker(A) \f$.
		 */
		template <typename RingLike>
		inline SparseMatrix<RingLike> SparseMatrixRREFKernel(
			SparseMatrix<RingLike>& A,
			Ring* R,
			SparsePivots& pivots,
			ComputeResources<RingLike>& options
		) {
			return SparseRREF::sparse_mat_rref_kernel<typename RingLike::dtype,INDEX>(A, pivots, R->ring, options.opt);
		};


		/**
		 * @brief Computes \f$ \ker(\RREF(A)) \f$.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param A Matrix.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * @param options Multithreaded computing options.
		 * 
		 * @returns \f$ \ker(\RREF(A)) \f$.
		 */
		template <typename RingLike>
		inline SparseMatrix<RingLike> SparseMatrixKernel(
			SparseMatrix<RingLike>& A,
			Ring* R,
			ComputeResources<RingLike>& options
		) {
			SparsePivots pivots = SparseMatrixRREF<RingLike>(A, R, options);
			return SparseRREF::sparse_mat_rref_kernel<typename RingLike::dtype,INDEX>(A, pivots, R->ring, options.opt);
		};
	};
}

#include "ATEAMS++/arithmetic/kernel.h"

#endif

