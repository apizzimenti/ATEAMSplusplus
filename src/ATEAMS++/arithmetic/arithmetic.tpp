
#ifndef ATEAMS_ARITHMETIC_T
#define ATEAMS_ARITHMETIC_T

#ifndef ATEAMS_ARITHMETIC_H
#error __FILE__ should only be included from arithmetic/arithmetic.h.
#endif

#include "ATEAMS++/arithmetic/arithmetic.h"

#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_vec.h>
#include <SparseRREF/scalar.h>

#include "ATEAMS++/arithmetic/options.h"

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
		INDEX maxindex = max(u(u.size()-1), v(v.size()-1));
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

		// Now, on each thread, copy each chunk of each sparse vector, then
		// add to `chunks`.
		for (int thread=0; thread < threads; thread++) {
			options.opt->pool.detach_task([thread, &u, &v, &R, &options] {
				int lo = options.parallel->indexBlocks[thread][0];
				int hi = options.parallel->indexBlocks[thread][1];
				SparseVector<RingLike> upart, vpart;

				// Get the relevant parts of u and v.
				for (int i=0; i < u.size(); i++) {
					if (lo <= u(i) && u(i) < hi) upart.push_back(u(i),u[i]);
				}

				for (int i=0; i < v.size(); i++) {
					if (lo <= v(i) && v(i) < hi) vpart.push_back(v(i),v[i]);
				}

				// Now, perform the addition, then write to the vector of chunks.
				sparse_vec_add<INDEX>(upart, vpart, R->ring);
				options.parallel->vectorBlocks[thread] = upart; // IS THIS WHERE THE COPY HAPPENS???
			});
		}
		options.opt->pool.wait();

		// Once we're done waiting, reconstitute the vector.
		u.zero();

		for (int thread=0; thread < threads; thread++) {
			for (int i=0; i < options.parallel->vectorBlocks[thread].size(); i++) {
				u.push_back(options.parallel->vectorBlocks[thread](i), options.parallel->vectorBlocks[thread][i]);
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
		if (options.parallel->enabled) return parallelSparseVectorAddition<Zp>(u, v, R, options);
		return serialSparseVectorAddition<Zp>(u, v, R);
	};


	template <>
	inline SparseVector<Z2> SparseVectorAddition<Z2>(
		SparseVector<Z2>& u,
		SparseVector<Z2>& v,
		Ring* R,
		arithmetic::ComputeOptions<Z2>& options
	) {
		if (options.parallel->enabled) return parallelSparseVectorAddition<Z2>(u, v, R, options);
		return serialSparseVectorAddition<Z2>(u, v, R);
	};
	
	
	// /**
	//  * @brief Sparse vector addition \f$ \vec u + \vec v \f$ for vectors over @ref Zp,
	//  * in parallel.
	//  * 
	//  * @param u Vector.
	//  * @param v Vector.
	//  * @param R (Pointer to) the coefficient ring @ref Zp.
	//  * @param options Parallel computing options.
	//  * 
	//  * @returns \f$ \vec w = \vec u + \vec v \f$.
	//  */
	// SparseVector<Zp> SparseVectorAddition(
	// 	SparseVector<Zp>& u,
	// 	SparseVector<Zp>& v,
	// 	Ring* R,
	// 	arithmetic::ComputeOptions<Zp>& options
	// ) {
	// 	// Get the number of threads, then the max/min of the largest/smallest
	// 	// indices in u and v.
	// 	int _threads = options.opt->pool.get_thread_count();
	// 	int threads = _threads*_threads;
	// 	INDEX maxindex = max(u(u.size()-1), v(v.size()-1));
	// 	INDEX minindex = min(u(0), v(0));

	// 	// Specify the index options.parallel->indexBlocks.
	// 	int mod = (maxindex-minindex)/threads;
	// 	int width = mod + 1;

	// 	vector<vector<int>> options.parallel->indexBlocks(threads, vector<int>(2,0));
	// 	options.parallel->indexBlocks[0][0] = minindex;
	// 	options.parallel->indexBlocks[0][1] = minindex+width;

	// 	for (int t=1; t < threads; t++) {
	// 		options.parallel->indexBlocks[t][0] = options.parallel->indexBlocks[t-1][1];
	// 		options.parallel->indexBlocks[t][1] = options.parallel->indexBlocks[t][0] + width;
	// 	}

	// 	options.parallel->indexBlocks[threads-1][1] = maxindex+1;

	// 	// Create a bucket for a vector on each thread; specify the length
	// 	// ahead-of-time, so that no vectors are overwriting one another.
	// 	vector<SparseVector<Zp>> chunks(threads);

	// 	// Now, on each thread, copy each chunk of each sparse vector, then
	// 	// add to `chunks`.
	// 	for (int thread=0; thread < threads; thread++) {
	// 		options.opt->pool.detach_task([thread, options.parallel->indexBlocks, &u, &v, &R, &chunks] {
	// 			SparseVector<Zp> upart, vpart;
	// 			int lo = options.parallel->indexBlocks[thread][0], hi = options.parallel->indexBlocks[thread][1];

	// 			// Get the relevant parts of u and v.
	// 			for (int i=0; i < u.size(); i++) {
	// 				if (lo <= u(i) && u(i) < hi) upart.push_back(u(i),u[i]);
	// 			}

	// 			for (int i=0; i < v.size(); i++) {
	// 				if (lo <= v(i) && v(i) < hi) vpart.push_back(v(i),v[i]);
	// 			}

	// 			// Now, perform the addition, then write to the vector of chunks.
	// 			sparse_vec_add<INDEX>(upart, vpart, R->ring);
	// 			chunks[thread] = upart;
	// 		});
	// 	}
	// 	options.opt->pool.wait();

	// 	// Once we're done waiting, reconstitute the vector.
	// 	u.zero();

	// 	for (int thread=0; thread < threads; thread++) {
	// 		for (int i=0; i < chunks[thread].size(); i++) {
	// 			u.push_back(chunks[thread](i), chunks[thread][i]);
	// 		}
	// 	}

	// 	return u;
	// };
}

#endif