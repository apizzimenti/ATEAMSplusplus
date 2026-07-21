
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
		 * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Zp.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) the coefficient ring @ref Zp.
		 * 
		 * @returns \f$ \vec w = \vec u + \vec v \f$.
		 */
		inline SparseVector<Zp> SparseVectorAddition(
			SparseVector<Zp>& u,
			SparseVector<Zp>& v,
			Ring* R
		) {
			sparse_vec_add<INDEX>(u, v, R->ring);
			return u;
		};
		
		
		/**
		 * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Zp,
		 * in parallel.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) the coefficient ring @ref Zp.
		 * @param options Parallel computing options.
		 * 
		 * @returns \f$ \vec w = \vec u + \vec v \f$.
		 */
		inline SparseVector<Zp> SparseVectorAddition(
			SparseVector<Zp>& u,
			SparseVector<Zp>& v,
			Ring* R,
			arithmetic::ComputeOptions& options
		) {
			// Get the number of threads, then the max/min of the largest/smallest
			// indices in u and v.
			int _threads = options.opt->pool.get_thread_count();
			int threads = _threads*_threads;
			INDEX maxindex = std::max(u(u.size()-1), v(v.size()-1));
			INDEX minindex = std::min(u(0), v(0));

			// Specify the index ranges.
			int mod = (maxindex-minindex)/threads;
			int width = mod + 1;

			std::vector<std::vector<int>> ranges(threads, std::vector<int>(2,0));
			ranges[0][0] = minindex;
			ranges[0][1] = minindex+width;

			for (int t=1; t < threads; t++) {
				ranges[t][0] = ranges[t-1][1];
				ranges[t][1] = ranges[t][0] + width;
			}

			ranges[threads-1][1] = maxindex+1;

			// Create a bucket for a vector on each thread; specify the length
			// ahead-of-time, so that no vectors are overwriting one another.
			std::vector<SparseVector<Zp>> chunks(threads);

			// Now, on each thread, copy each chunk of each sparse vector, then
			// add to `chunks`.
			for (int thread=0; thread < threads; thread++) {
				options.opt->pool.detach_task([thread, ranges, &u, &v, &R, &chunks] {
					SparseVector<Zp> upart, vpart;
					int lo = ranges[thread][0], hi = ranges[thread][1];

					// Get the relevant parts of u and v.
					for (int i=0; i < u.size(); i++) {
						if (lo <= u(i) && u(i) < hi) upart.push_back(u(i),u[i]);
					}

					for (int i=0; i < v.size(); i++) {
						if (lo <= v(i) && v(i) < hi) vpart.push_back(v(i),v[i]);
					}

					// Now, perform the addition, then write to the vector of chunks.
					sparse_vec_add<INDEX>(upart, vpart, R->ring);
					chunks[thread] = upart;
				});
			}
			options.opt->pool.wait();

			// Once we're done waiting, reconstitute the vector.
			u.zero();

			for (int thread=0; thread < threads; thread++) {
				for (int i=0; i < chunks[thread].size(); i++) {
					u.push_back(chunks[thread](i), chunks[thread][i]);
				}
			}

			return u;
		};

		// /**
		//  * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Z2.
		//  * 
		//  * @param u Vector.
		//  * @param v Vector.
		//  * @param R (Pointer to) the coefficient ring @ref Z2.
		//  * 
		//  * @returns \f$ \vec w = \vec u + \vec v \f$.
		//  */
		// inline SparseVector<Z2> SparseVectorAddition(
		// 	SparseVector<Z2>& u,
		// 	SparseVector<Z2> v,
		// 	Ring* R
		// ) {
		// 	sparse_vec_add<INDEX>(u, v, R->ring);
		// 	return u;
		// };

		/**
		 * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Q.
		 * 
		 * @param u Vector.
		 * @param v Vector.
		 * @param R (Pointer to) the coefficient ring @ref Q.
		 * 
		 * @returns \f$ \vec w = \vec u + \vec v \f$.
		 */
		inline SparseVector<Q> SparseVectorAddition(
			SparseVector<Q>& u,
			SparseVector<Q>& v,
			Ring* R,
			arithmetic::ComputeOptions& options
		) {
			sfmpq_vec_addsub_mul<INDEX,false>(u, v, (Q::dtype)1);
			return u;
		};


		/**
		 * @brief Sparse vector rescaling \f$ a \vec x \f$.
		 * @tparam RingLike A coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @param a Scalar.
		 * @param x Vector.
		 * @param R (Pointer to) a coefficient @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @returns \f$ a \vec x \f$.
		 */
		template <typename RingLike>
		inline SparseVector<RingLike> SparseVectorRescaling(
			typename RingLike::dtype a,
			SparseVector<RingLike>& x,
			Ring* R
		) {
			sparse_vec_rescale<INDEX,typename RingLike::dtype>(x, a, R->ring);
			return x;
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
			ComputeOptions& options
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
			ComputeOptions& options
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
			ComputeOptions& options
		) {
			SparsePivots pivots = SparseMatrixRREF<RingLike>(A, R, options);
			return SparseRREF::sparse_mat_rref_kernel<typename RingLike::dtype,INDEX>(A, pivots, R->ring, options.opt);
		};
	};
}

#include "ATEAMS++/arithmetic/kernel.h"

#endif

