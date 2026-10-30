
#ifndef ATEAMS_ARITHMETIC_H
#define ATEAMS_ARITHMETIC_H

#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_vec.h>
#include <SparseRREF/scalar.h>

#include "ATEAMS++/arithmetic/resources.h"

namespace ATEAMS {
	/**
	 * @brief Sparse, finite-field linear algebra and support routines.
	 */
	namespace arithmetic {

		/**
		 * @brief Sparse vector dot product \f$u \cdot v\f$ for vectors over a
		 * @ref ATEAMS::Ring.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) a coefficient @ref Ring \f$R\f$, like @ref Zp or @ref Q.
		 * 
		 * @returns \f$(u \cdot v) \in R\f$.
		 */
		template <typename RingLike>
		inline typename RingLike::dtype SparseVectorDot(
			SparseVector<RingLike>& u,
			const SparseVector<RingLike>& v,
			Ring* R
		) {
			return sparse_vec_dot<INDEX, typename RingLike::dtype>(u, v, R->ring);
		};

		/**
		 * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over a
		 * @ref ATEAMS::Ring
	 	 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * @param resources Computing resources with cached arithmetic tables.
		 */
		template <typename RingLike>
		inline void SparseVectorAddition(
			SparseVector<RingLike>& u,
			SparseVector<RingLike>& v,
			Ring* R,
			ComputeResources<RingLike>& resources
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
					auto entry = resources.arithmetic->add[u[left-1]][v[right-1]];

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
					u[out-1] = v[right-1];

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
					u[out-1] = u[left-1];
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
			#pragma omp parallel for default(shared) schedule(static)
			for (size_t i=left; i < out; i++) u[i] = 0;

			// Canonicalize.
			u.canonicalize();
		}

		/**
		 * @brief Template specialization; sparse vector addition \f$ \vec u + \vec v \f$
		 * for vectors over @ref ATEAMS::Z2.
	 	 * @tparam None.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * @param resources Computing resources with cached arithmetic tables.
		 */
		template <>
		inline void SparseVectorAddition<Z2>(
			SparseVector<Z2>& u,
			SparseVector<Z2>& v,
			Ring* R,
			ComputeResources<Z2>& resources
		) {
			sparse_vec_add<INDEX>(u, v, R->ring);
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
			return sparse_mat_dot_sparse_vec<typename RingLike::dtype,INDEX>(
				A, x, R->ring
			);
		};

		// /**
		//  * @brief Sparse (right) matrix multiplication \f$ A \vec x \f$, performed
		//  * 	by computing dot products \f$A_i \cdot \vec x\f$ in parallel. Uses
		//  *  OpenMP.
		//  * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		//  * 
		//  * @param A Matrix.
		//  * @param x Vector.
		//  * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		//  * @param resources Compute resources.
		//  * 
		//  * @returns \f$ \vec y = A \vec x \f$.
		//  */
		// template <typename RingLike>
		// inline SparseVector<RingLike> SparseRightMultiplication(
		// 	SparseMatrix<RingLike>& A,
		// 	SparseVector<RingLike>& x,
		// 	Ring* R,
		// 	ComputeResources<RingLike>& resources
		// ) {
		// 	// If parallelism isn't enabled, kick to the other version.
		// 	if (!resources.parallel->enabled) return SparseRightMultiplication<RingLike>(A, x, R);

		// 	std::vector<SparseVector<RingLike>> chunks(omp_get_max_threads(), SparseVector<RingLike>());

		// 	// Compute the dot products in parallel.
		// 	#pragma omp parallel for default(shared) schedule(static) if (A.rows.size() > 1e5)
		// 	for (int i=0; i < A.rows.size(); i++) {
		// 		int tid = omp_get_thread_num();
		// 		typename RingLike::dtype r = arithmetic::SparseVectorDot<RingLike>(A.rows[i], x, R);

		// 		if (r > 0) resources.typed->parts[tid].push_back((INDEX)i, r);
		// 	}

		// 	// Reconstitute the vector.
		// 	SparseVector<RingLike> r;

		// 	for (int j=0; j < resources.typed->parts.size(); j++) {
		// 		// Skip this row if there's nothing in it.
		// 		if (resources.typed->parts[j].size() < 1) continue;

		// 		for (int k=0; k < resources.typed->parts[j].size(); k++) {
		// 			r.push_back(resources.typed->parts[j](k), resources.typed->parts[j][k]);
		// 		}
		// 	}
		// 	r.sort_indices();
		// 	resources.typed->flush();

		// 	return r;
		// }


		/**
		 * @brief Sparse (right) matrix multiplication \f$ A \vec x \f$, performed
		 * 	by computing dot products \f$A_i \cdot \vec x\f$ in parallel. Uses
		 *  `BS::thread_pool`.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param A Matrix.
		 * @param x Vector.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * @param resources Compute resources.
		 * 
		 * @returns \f$ \vec y = A \vec x \f$.
		 */
		template <typename RingLike>
		inline SparseVector<RingLike> SparseRightMultiplication(
			SparseMatrix<RingLike>& A,
			SparseVector<RingLike>& x,
			Ring* R,
			ComputeResources<RingLike>& resources
		) {
			// If parallelism isn't enabled, kick to the other version.
			if (!resources.parallel->enabled) return SparseRightMultiplication<RingLike>(A, x, R);

			// If we haven't yet decided on a preferred number of blocks, let's
			// compute that now.
			if (resources.typed->blocks < 0) {
				resources.typed->defaultPreferredBlocks(A.rows.size());
				resources.typed->build();
			}

			// Compute the dot products in parallel.
			resources.opt->pool.detach_blocks(
				0, A.rows.size(),
				[&](const size_t start, const size_t end) {
					if (resources.opt->abort) return;

					// Get the thread ID
					auto tid = SparseRREF::thread_id();
					typename RingLike::dtype r;

					// Loop over each block.
					for (size_t i=start; i < end; i++) {
						r = arithmetic::SparseVectorDot<RingLike>(A.rows[i], x, R);
						if (r > 0) resources.typed->parts[tid].push_back((INDEX)i, r);
					}
				},
				resources.typed->blocks
			);
			resources.opt->pool.wait();

			// Reconstitute the vector.
			SparseVector<RingLike> r;

			for (int j=0; j < resources.typed->parts.size(); j++) {
				// Skip this row if there's nothing in it.
				if (resources.typed->parts[j].size() < 1) continue;

				for (int k=0; k < resources.typed->parts[j].size(); k++) {
					r.push_back(resources.typed->parts[j](k), resources.typed->parts[j][k]);
				}
			}
			r.sort_indices();
			resources.typed->flush();

			return r;
		}

		/**
		 * @brief Given a (sparse) vector \f$\vec x\f$ and scalar \f$a\f$, computes
		 * 	\f$a \vec x\f$.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param a Scalar..
		 * @param x Vector.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * @param resources Compute resources.
		 */
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
		 * @param resources Multithreaded computing resources.
		 * 
		 * @returns \f$ \RREF(A) \f$.
		 */
		template <typename RingLike>
		inline SparsePivots SparseMatrixRREF(
			SparseMatrix<RingLike>& A,
			Ring* R,
			ComputeResources<RingLike>& resources
		) {
			return SparseRREF::sparse_mat_rref<typename RingLike::dtype,INDEX>(A, R->ring, resources.opt);
		};

		/**
		 * @brief Computes \f$ \ker(A) \f$, where \f$A\f$ is in reduced row
		 * 	echelon form.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param A Matrix.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * @param pivots Pivot indices of \f$ A \f$.
		 * @param resources Multithreaded computing resources.
		 * 
		 * @returns \f$ \ker(A) \f$.
		 */
		template <typename RingLike>
		inline SparseMatrix<RingLike> SparseMatrixRREFKernel(
			SparseMatrix<RingLike>& A,
			Ring* R,
			SparsePivots& pivots,
			ComputeResources<RingLike>& resources
		) {
			return SparseRREF::sparse_mat_rref_kernel<typename RingLike::dtype,INDEX>(A, pivots, R->ring, resources.opt);
		};


		/**
		 * @brief Computes \f$ \ker(\RREF(A)) \f$.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param A Matrix.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * @param resources Multithreaded computing resources.
		 * 
		 * @returns \f$ \ker(\RREF(A)) \f$.
		 */
		template <typename RingLike>
		inline SparseMatrix<RingLike> SparseMatrixKernel(
			SparseMatrix<RingLike>& A,
			Ring* R,
			ComputeResources<RingLike>& resources
		) {
			SparsePivots pivots = SparseMatrixRREF<RingLike>(A, R, resources);
			return SparseRREF::sparse_mat_rref_kernel<typename RingLike::dtype,INDEX>(A, pivots, R->ring, resources.opt);
		};
	};
}

#include "ATEAMS++/arithmetic/kernel.h"

#endif

