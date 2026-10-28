
#ifndef ATEAMS_MODELS_SWENDSENWANG_T
#define ATEAMS_MODELS_SWENDSENWANG_T

#ifndef ATEAMS_MODELS_SWENDSENWANG_H
#error __FILE__ should only be included from models/SwendsenWang.h.
#endif

#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/SwendsenWang.h"
#include "ATEAMS++/arithmetic/arithmetic.h"

#include <SparseRREF/sparse_vec.h>

#include <random>
#include <algorithm>
#include <cassert>
#include <execution>


namespace ATEAMS {
	namespace models {

		template <typename RingLike>
		SwendsenWang<RingLike>::State SwendsenWang<RingLike>::sample(
			int t,
			SwendsenWang<RingLike>::State& state,
			arithmetic::ComputeResources<RingLike>& resources
		) {
			// Compute the temperature and probability of including particular (d-1)-cells.
			double temp = this->temperatureFunction(t),
				p = 1-exp(temp);

			// If we're debugging, check that the inclusion probability is actually within
			// acceptable parameters.
			if (this->DEBUG) {
				cerr << std::format("verifying 0 ≤ {} ≤ 1", p) << endl;
				assert((0 <= p) && (p <= 1));
			}

			int d = this->dimension,
				N = this->complex->Cells[d],
				satisfied;

			// Perform independent Bernoulli trials over each (d-1)-cell to decide whether
			// it gets included in the computation. We do this in steps: given a
			// cochain f and the dth coboundary matrix B,
			//
			//		1.	perform the multiplication Bf (by performing the dot
			//			products B(k)•f, where B(k) is the kth row of B, in
			//			parallel) and check which entries of Bf are nonzero;
			//		2.	iterate over all cells and flip a density-p coin (i.e.
			//			a coin that shows "heads" with probability p) and exclude
			//			each cell for which the coin shows "tails." (For any cell
			//			already excluded, this does nothing; for any cell not yet
			//			excluded, this excludes it.)
			//
			// Passing `resources` as a fourth argument performs the underlying
			// matrix multiplication in parallel.

			vector<int> exclude = arithmetic::nonzeroIndices<RingLike>(
				this->complex->Coboundary.Matrices[d],
				state.cochain,
				this->coefficients,
				resources
			);

			satisfied = N-exclude.size();

			// If we're debugging, verify that the cells corresponding to indices
			// reported here are actually not sent to 0.
			if (this->DEBUG) {
				for (auto& e : exclude) {
					typename RingLike::dtype r = arithmetic::SparseVectorDot<RingLike>(
						this->complex->Coboundary.Matrices[d][e],
						state.cochain,
						this->coefficients
					);
					assert(r > 0);
				}
			}

			// We then iterate over all cells and flip a density-p coin (i.e.
			// a coin that shows "heads" with probability p) and exclude each cell
			// for which the coin shows "tails." (For any cell already excluded,
			// this does nothing; for any cell not yet excluded, this excludes it.)
			for (int i=0; i<N; i++) {
				if (this->unituniform(this->RNG) > p) exclude.push_back(i);
			}

			// Delete duplicates.
			if (N > 1e6) std::sort(std::execution::par, exclude.begin(), exclude.end());
			else std::sort(exclude.begin(), exclude.end());
			
			exclude.erase(std::unique(exclude.begin(), exclude.end()), exclude.end());

			// Sample from the kernel.
			SparseVector<RingLike> sample = arithmetic::submatrixKernelSample<RingLike>(
				this->complex->Coboundary.Matrices[d],	// complete dth coboundary matrix
				this->coefficients,						// ring
				exclude,								// rows to exclude
				this->intuniform,						// uniform random over field
				this->RNG,								// RNG
				resources,								// compute resources
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
				// there shouldn't be (i.e. it should have size 0). Adding a fourth
				// argument (passing `resources`) performs this multiplication in
				// parallel.
				SparseVector<RingLike> outcome = arithmetic::SparseRightMultiplication<RingLike>(
					cbd,
					sample,
					this->coefficients,
					resources
				);

				assert(outcome.size() < 1);
			}

			// Figure out which indices we're *in*cluding. (Is this necessary?)
			std::vector<int> included;
			
			set_difference(
				this->include.begin(), this->include.end(),
				exclude.begin(), exclude.end(),
				inserter(included, included.begin())
			);

			state.cochain = sample;
			state.includes = included;
			state.satisfied = satisfied;				// TODO satisfied of the previous cochain
			state.occupied = N-exclude.size();
			state.energy = -exclude.size();				// TODO energy of the previous cochain
			state.t = t;

			return state;
		}

		template <typename RingLike>
		SwendsenWang<RingLike>::State SwendsenWang<RingLike>::initialize(
			SwendsenWang<RingLike>::State& state
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
		SwendsenWang<RingLike>::State SwendsenWang<RingLike>::initialize(
			SparseVector<RingLike> c,
			SwendsenWang<RingLike>::State& state
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
			std::vector<int> include(this->complex->Cells[this->dimension]);
			std::iota(include.begin(), include.end(), 0);
			this->include = include;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->unituniform = std::uniform_real_distribution<double>(0,1);

			this->intuniform = std::uniform_int_distribution<int>(
				0,
				this->coefficients->characteristic > 0 ? this->coefficients->characteristic-1 : 1
			);
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
			std::vector<int> include(this->complex->Cells[this->dimension]);
			std::iota(include.begin(), include.end(), 0);
			this->include = include;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->unituniform = std::uniform_real_distribution<double>(0,1);

			this->intuniform = std::uniform_int_distribution<int>(
				0,
				this->coefficients->characteristic > 0 ? this->coefficients->characteristic-1 : 1
			);
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
			std::vector<int> include(this->complex->Cells[this->dimension]);
			std::iota(include.begin(), include.end(), 0);
			this->include = include;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->unituniform = std::uniform_real_distribution<double>(0,1);

			this->intuniform = std::uniform_int_distribution<int>(
				0,
				this->coefficients->characteristic > 0 ? this->coefficients->characteristic-1 : 1
			);
		}
	}
}


#endif
