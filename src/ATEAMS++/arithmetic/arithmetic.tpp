
#ifndef ATEAMS_ARITHMETIC_T
#define ATEAMS_ARITHMETIC_T

#ifndef ATEAMS_ARITHMETIC_H
#error __FILE__ should only be included from arithmetic/arithmetic.h.
#endif

#include "ATEAMS++/arithmetic/arithmetic.h"
#include "ATEAMS++/arithmetic/options.h"

#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_vec.h>
#include <SparseRREF/scalar.h>

using namespace SparseRREF;
using namespace std;

namespace ATEAMS::arithmetic {

	template <typename RingLike>
	inline SparseVector<RingLike> serialSparseVectorAddition(
		SparseVector<RingLike>& u,
		SparseVector<RingLike>& v,
		Ring* R
	) {
		sparse_vec_add<INDEX>(u, v, R->ring);
		return u;
	};

	template <typename RingLike>
	inline SparseVector<RingLike> parallelSparseVectorAddition(
		SparseVector<RingLike>& u,
		SparseVector<RingLike>& v,
		Ring* R,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		// Get the number of threads, then the max/min of the largest/smallest
		// indices in u and v.
		int threads = options.parallel->threads;
		INDEX maxindex = max(u(u._nnz-1), v(v._nnz-1));
		INDEX minindex = min(u(0), v(0));

		// Specify index blocks for this run.
		int mod = (maxindex-minindex)/threads;
		int width = mod + 1;

		options.parallel->indexBlocks[0][0] = minindex;
		options.parallel->indexBlocks[0][1] = minindex+width;

		for (int t=1; t < threads; t++) {
			options.parallel->indexBlocks[t][0] = options.parallel->indexBlocks[t-1][1];
			options.parallel->indexBlocks[t][1] = options.parallel->indexBlocks[t][0] + width;
		}

		options.parallel->indexBlocks[threads-1][1] = maxindex+1;

		// I think parallelizing in this way creates a ton of overhead just from
		// submitting to the thread pool. May need to re-think.
		for (int thread=0; thread < threads; thread++) {
			options.opt->pool.detach_task([thread, &u, &v, &R, &options] {
				// Determine which blocks of indices we're working with.
				int lo = options.parallel->indexBlocks[thread][0];
				int hi = options.parallel->indexBlocks[thread][1];

				// Follow a
				//
				//	1. reserve max space,
				//	2. overwrite entries,
				//	4. resize
				//
				// pattern for each vector in the scratch space. Hopefully this
				// saves us some memory, since we'll be re-sizing fewer times
				// than if we used push_back. Not calling .reserve() a second time
				// (to cut the vector down to size) may have some effects, but
				// we'll see what the profiling data says.
				options.parallel->lScratch[thread].reserve(hi-lo, false);
				options.parallel->rScratch[thread].reserve(hi-lo, false);
				int lnnz = 0, rnnz = 0;

				for (int i=0; i < max(u._nnz, v._nnz); i++) {
					// Check for entries in the first vector.
					if (i < u._nnz && lo <= u(i) && u(i) < hi) {
						options.parallel->lScratch[thread].indices[lnnz] = u(i);
						options.parallel->lScratch[thread].entries[lnnz] = u[i];
						lnnz++;
					}

					// Check for entries in the second vector.
					if (i < v._nnz && lo <= v(i) && v(i) < hi) {
						options.parallel->rScratch[thread].indices[rnnz] = v(i);
						options.parallel->rScratch[thread].entries[rnnz] = v[i];
						rnnz++;
					}
				}

				// Resize vectors.
				options.parallel->lScratch[thread]._nnz = lnnz;
				options.parallel->rScratch[thread]._nnz = rnnz;

				// Add.
				sparse_vec_add(
					options.parallel->lScratch[thread],
					options.parallel->rScratch[thread],
					R->ring
				);
			});
		}
		options.opt->pool.wait();

		// Since the number of threads is constant and the call to .size() is
		// constant, this is constant-time. Maybe use ._nnz instead?
		int nnz = 0;
		for (int thread=0; thread < threads; thread++) {
			nnz += options.parallel->lScratch[thread]._nnz;
		}

		// Zero the vector, then re-size it so we can write directly into it
		// rather than pushing back.
		u.reserve(nnz, false);
		u._nnz = nnz;

		int t = 0;

		for (int thread=0; thread < threads; thread++) {
			for (int i=0; i < options.parallel->lScratch[thread]._nnz; i++) {
				u.indices[t] = options.parallel->lScratch[thread](i);
				u.entries[t] = options.parallel->lScratch[thread][i];
				t++;
			}
		}

		return u;
	};


	template <>
	inline SparseVector<Zp> SparseVectorAddition<Zp>(
		SparseVector<Zp>& u,
		SparseVector<Zp>& v,
		Ring* R,
		arithmetic::ComputeOptions<Zp>& options
	) {
		// if (options.parallel->enabled) return parallelSparseVectorAddition<Zp>(u, v, R, options);
		// return serialSparseVectorAddition<Zp>(u, v, R);
		return serialSparseVectorAddition<Zp>(u, v, R);
	};


	template <>
	inline SparseVector<Z2> SparseVectorAddition<Z2>(
		SparseVector<Z2>& u,
		SparseVector<Z2>& v,
		Ring* R,
		arithmetic::ComputeOptions<Z2>& options
	) {
		// if (options.parallel->enabled) return parallelSparseVectorAddition<Z2>(u, v, R, options);
		// return serialSparseVectorAddition<Z2>(u, v, R);
		return serialSparseVectorAddition<Z2>(u, v, R);
	};
}

#endif