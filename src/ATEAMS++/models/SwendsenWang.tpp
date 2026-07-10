
#ifndef ATEAMS_MODELS_SWENDSENWANG_T
#define ATEAMS_MODELS_SWENDSENWANG_T

#ifndef ATEAMS_MODELS_SWENDSENWANG_H
#error __FILE__ should only be included from models/SwendsenWang.h.
#endif

#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/SwendsenWang.h"
#include "ATEAMS++/arithmetic/kernel.h"

#include <SparseRREF/sparse_mat.h>

#include <random>
#include <algorithm>
#include <cassert>


namespace ATEAMS {
	namespace models {

		template <typename RingLike>
		ModelState<RingLike,SparseVector> SwendsenWang<RingLike>::sample(
			int t,
			ModelState<RingLike,SparseVector>& state,
			arithmetic::ThreadOptions& options
		) {
			// Compute the temperature and probability of including particular (d-1)-cells.
			double temp = this->temperatureFunction(t);
			double p = 1-exp(temp);

			// If we're debugging, check that the inclusion probability is actually within
			// acceptable parameters.
			if (this->DEBUG) {
				cerr << std::format("verifying 0 ≤ {} ≤ 1", p) << endl;
				assert((0 <= p) && (p <= 1));
			}

			// Perform independent Bernoulli trials over each (d-1)-cell to decide whether
			// it gets included in the computation. We do this by performing a (sparse!) matrix
			// multiplication, then checking which of the entries are zero.
			size_t d = this->dimension;

			SparseVector<RingLike> coeffs = sparse_mat_dot_sparse_vec<typename RingLike::dtype,INDEX>(
				this->complex->Coboundary.Matrices[d],
				state.cochain,
				this->coefficients->ring
			);
			coeffs.compress();

			// The zero entries are those that don't show up in the indices, so we
			// automatically exclude these.
			int N = this->complex->Cells[d];
			std::set<size_t> exclude;
			
			if (coeffs.size() > 0) {
				for (int e=0; e < coeffs.size(); e++) exclude.insert((size_t)e);
			}

			// We then iterate over the remaining indices and flip a weighted coin (i.e.
			// perform a density-p Bernoulli trial) to decide whether to keep that
			// particular row (corresponding to a d-cell on which the current cochain
			// evaluates to 0).
			for (int i=0; i<N; i++) {
				if (this->unituniform(this->RNG) > p) exclude.insert(i);
			}

			// Figure out which indices we're *in*cluding.
			std::set<int> includeHere;
			
			set_difference(
				this->include.begin(), this->include.end(),
				exclude.begin(), exclude.end(),
				inserter(includeHere, includeHere.begin())
			);

			// Now, sample from the kernel.
			SparseVector<RingLike> sample = arithmetic::submatrixKernelSample<RingLike>(
				this->complex->Coboundary.Matrices[d],	// complete dth coboundary matrix
				this->coefficients,								// ring
				exclude,								// rows to exclude
				this->intuniform,						// uniform random over field
				this->RNG,								// RNG
				options,								// compute options
				this->DEBUG								// debugging
			);

			// If we're debugging, check whether we actually sampled from the kernel.
			if (this->DEBUG) {
				// Copy and pare down the coboundary matrix.
				SparseMatrix<RingLike> cbd = this->complex->Coboundary.Matrices[d];
				for (auto i : exclude) cbd[i].zero();
				cbd.clear_zero_row();
				cbd.compress();

				// Multiply, and check whether there's anything in the resulting vector;
				// there shouldn't be (i.e. it should have size 0).
				SparseVector<RingLike> outcome = sparse_mat_dot_sparse_vec<typename RingLike::dtype,INDEX>(
					cbd,
					sample,
					this->coefficients->ring
				);

				assert(outcome.size() < 1);
			}

			state.cochain = sample;
			state.includes = std::vector<int>(includeHere.begin(), includeHere.end());
			state.t = t;

			return state;
		}

		template <typename RingLike>
		ModelState<RingLike,SparseVector> SwendsenWang<RingLike>::initialize(
			ModelState<RingLike,SparseVector>& state
		) {
			size_t dimension = this->dimension-1;
			int N = this->complex->Cells[dimension];

			SparseVector<RingLike> cochain;
			for (int i=0; i < N; i++) cochain.push_back(
				(INDEX)i, (typename RingLike::dtype)this->intuniform(this->RNG)
			);

			cochain.compress();
			state.cochain = cochain;
			return state;
		}


		template <typename RingLike>
		ModelState<RingLike,SparseVector> SwendsenWang<RingLike>::initialize(
			SparseVector<RingLike> c,
			ModelState<RingLike,SparseVector>& state
		) {
			state.cochain = c;
			return state;
		}


		template <typename RingLike>
		SwendsenWang<RingLike>::SwendsenWang(
			complexes::Complex<RingLike>* complex,
			Ring* R,
			int dimension,
			bool DEBUG
		) : Model<RingLike,SparseVector>(R, dimension, DEBUG) {
			this->complex = complex;
			this->temperatureFunction = statistics::selfdual(this->coefficients);

			// Determine the field and build the boundary matrices for the Complex.
			this->complex->constructBoundaryMatrices(this->coefficients);

			// Default set of indices to compare against.
			std::set<int> _include;
			for (int t=0; t < this->complex->Cells[this->dimension]; t++) _include.insert(t);

			this->include = _include;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->unituniform = std::uniform_real_distribution<double>(0,1);

			int mod = (int)R->characteristic;
			this->intuniform = std::uniform_int_distribution<int>(0, mod > 0 ? mod : 1);
		};

		template <typename RingLike>
		SwendsenWang<RingLike>::SwendsenWang(
			complexes::Complex<RingLike>* complex,
			Ring* R,
			int dimension,
			std::function<double(int)> temperatureFunction,
			bool DEBUG
		) : Model<RingLike,SparseVector>(R, dimension, DEBUG) {
			this->complex = complex;
			this->temperatureFunction = temperatureFunction;

			// Determine the field and build the boundary matrices for the Complex.
			this->complex->constructBoundaryMatrices(this->coefficients);

			// Default set of indices to compare against.
			std::set<int> _include;
			for (int t=0; t < this->complex->Cells[this->dimension]; t++) _include.insert(t);

			this->include = _include;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->unituniform = std::uniform_real_distribution<double>(0,1);

			int mod = (int)R->characteristic;
			this->intuniform = std::uniform_int_distribution<int>(0, mod > 0 ? mod : 1);
		};

		template <typename RingLike>
		SwendsenWang<RingLike>::SwendsenWang(
			complexes::Complex<RingLike>* complex,
			ModelParameters parameters
		) : Model<RingLike,SparseVector>(
			parameters.coefficients,
			parameters.dimension,
			parameters.DEBUG
		) {
			this->complex = complex;
			this->temperatureFunction = parameters.temperatureFunction;

			// Determine the field and build the boundary matrices for the Complex.
			this->complex->constructBoundaryMatrices(this->coefficients);

			// Default set of indices to compare against.
			std::set<int> _include;
			for (int t=0; t < this->complex->Cells[this->dimension]; t++) _include.insert(t);

			this->include = _include;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->unituniform = std::uniform_real_distribution<double>(0,1);

			int mod = (int)this->coefficients->characteristic;
			this->intuniform = std::uniform_int_distribution<int>(0, mod > 0 ? mod : 1);
		}
	}
}


#endif
