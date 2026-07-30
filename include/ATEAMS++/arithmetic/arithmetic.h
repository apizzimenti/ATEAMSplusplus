
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
			const SparseVector<RingLike>& v,
			Ring* R,
			ComputeResources<RingLike>& options
		) {
			// If one of the vectors is zero, just return.
			if (v.nnz() == 0) return;
			if (u.nnz() == 0) return;

			// Start from the largest index of each vector.
			size_t left = u.nnz();
			size_t right = v.nnz();
			size_t out = u.nnz() + v.nnz();

			// Preallocate enough space to accommodate all entries in both vectors.
			if (u.alloc() < out) u.reserve(out);
			u.resize(out);

			while (left > 0 && right > 0) {
				// If we're pointing to the same index...
				if (u(left-1) == v(right-1)) {
					auto entry = options.arithmetic->add[(int)u[left-1]][(int)v[right-1]];

					// ... and we get a nonzero entry, then the index of this entry
					// is the same as in the original vector, and the entry is the
					// result of the arithmetic; we then decrement the `out` pointer...
					if (entry != 0) {
						u(out-1) = u(left-1);
						u[out-1] = entry;
						out--;
					}

					// ... and decrement both left and right pointers, since we
					// found a result in both.
					left--;
					right--;
				}
				
				// If we *aren't* pointing at the same index, and the index on the
				// left is smaller than the index on the right...
				else if (u(left-1) < v(right-1)) {
					// ... then we set the index of this entry in the output vector
					// to the index of the one on the right (since this preserves
					// the order), then set the entry...
					u(out-1) = v(right-1);
					u[out-1] = (typename RingLike::dtype)v[right-1];

					// ... then decrement the right pointer and the out pointer,
					// since they track together in this case.
					right--;
					out--;
				}

				// Otherwise, the index on the right is smaller than the one on
				// the left, and we do the same thing as the "else if" from above
				// but with the roles reversed.
				else {
					u(out-1) = u(left-1);
					u[out-1] = (typename RingLike::dtype)u[left-1];
					left--;
					out--;
				}
			}

			// Then, while the right pointer is still bigger than zero, write
			// the remaining entries into the vector.
			while (right > 0) {
				u(out-1) = v(right-1);
				u[out-1] = v[right-1];
				right--;
				out--;
			}

			// Zero out the remaining entries, if there are any.
			for (size_t i=left; i < out; i++) u[i] = 0;

			// Canonicalize.
			u.canonicalize();
		}

		// /**
		//  * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Z2.
		//  * 
		//  * @param u Vector.
		//  * @param v Vector.
		//  * @param R (Pointer to) the coefficient ring @ref Q.
		//  */
		// template <>
		// inline void SparseVectorAddition<Z2>(
		// 	SparseVector<Z2>& u,
		// 	SparseVector<Z2>& v,
		// 	Ring* R,
		// 	ComputeResources<Z2>& options
		// ) {
		// 	SparseVectorAddition<Zp>(u, v, R, options);
		// };

		// /**
		//  * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Zp.
		//  * 
		//  * @param u Vector.
		//  * @param v Vector.
		//  * @param R (Pointer to) the coefficient ring @ref Q.
		//  */
		// template <>
		// inline void SparseVectorAddition<Zp>(
		// 	SparseVector<Zp>& u,
		// 	SparseVector<Zp>& v,
		// 	Ring* R,
		// 	ComputeResources<Zp>& options
		// ) {
			// sparse_vec_add<INDEX>(u, v, R->ring);
		// };

		// /**
		//  * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Q.
		//  * 
		//  * @param u Vector.
		//  * @param v Vector.
		//  * @param R (Pointer to) the coefficient ring @ref Q.
		//  */
		// template <>
		// inline void SparseVectorAddition<Q>(
		// 	SparseVector<Q>& u,
		// 	SparseVector<Q>& v,
		// 	Ring* R,
		// 	ComputeResources<Q>& options
		// ) {
		// 	sfmpq_vec_addsub_mul<INDEX,false>(u, v, (Q::dtype)1);
		// }

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


		// /**
		//  * @brief Sparse vector rescaling \f$ a \vec x \f$.
		//  * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		//  * 
		//  * @param a Scalar.
		//  * @param x Vector.
		//  * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		//  */
		// template <typename RingLike>
		// inline void SparseVectorRescaling(
		// 	typename RingLike::dtype a,
		// 	SparseVector<RingLike>& x,
		// 	Ring* R
		// ) {
		// 	sparse_vec_rescale<INDEX,typename RingLike::dtype>(x, a, R->ring);
		// };

		template <typename RingLike>
		inline void SparseVectorRescaling(
			typename RingLike::dtype a,
			SparseVector<RingLike>& x,
			Ring* R,
			ComputeResources<RingLike>& resources
		) {
			for (int t=0; t < x.size(); t++) {
				x[t] = resources.arithmetic->multiply[x[t]][a];
			}
		}

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

