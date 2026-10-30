
#ifndef ATEAMS_MODELS_INVADEDCLUSTER_T
#define ATEAMS_MODELS_INVADEDCLUSTER_T

#ifndef ATEAMS_MODELS_INVADEDCLUSTER_H
#error __FILE__ should only be included from models/InvadedCluster.h.
#endif

#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/InvadedCluster.h"
#include "ATEAMS++/arithmetic/kernel.h"
#include "ATEAMS++/topology/persistence.h"

#include <SparseRREF/sparse_vec.h>

#include <random>
#include <algorithm>
#include <cassert>


namespace ATEAMS {
	namespace models {

		template <typename RingLike>
		InvadedCluster<RingLike>::State InvadedCluster<RingLike>::sample(
			int t,
			InvadedCluster<RingLike>::State& state,
			arithmetic::ComputeResources<RingLike>& resources
		) {
			int d = this->dimension,
				N = this->complex->Cells[d];

			// These are the cells on which the cochain does not vanish.
			vector<int> exclude = arithmetic::nonzeroIndices<RingLike>(
				this->complex->Coboundary.Matrices[d],
				state.cochain,
				this->coefficients
			);

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

			// The cells to include are the ones not excluded.
			vector<int> include;

			set_difference(
				this->include.begin(), this->include.end(),
				exclude.begin(), exclude.end(),
				inserter(include, include.begin())
			);

			std::shuffle(include.begin(), include.end(), this->RNG);

			// Construct the filtration.
			std::vector<size_t> inside(include.begin(), include.end());
			inside.insert(inside.end(), exclude.begin(), exclude.end());

			int stop = this->complex->Breaks[d][1];
			int start = this->complex->Breaks[d][0];
			int offset = (d > 0) ? this->complex->Offsets[d-1] : 0;
			
			// Fill in all the indices for cells of dimension not equal to d.
			// TODO a small optimization: we only have to do the persistence algorithm
			// through the (d+1)-dimensional cells to capture the d-dimensional events.
			// This will probably save a boatload of time!
			for (int j=0; j < inside.size(); j++) this->filtration[j+start] = inside[j]+start;

			// Persist.
			std::vector<int> essential = topology::persistence::persistence<RingLike>(
				this->complex,
				this->filtration,
				this->coefficients,
				this->dimension,
				resources
			);

			if (this->DEBUG) {
				cerr << "filtration:" << endl;
				printvector<int>(this->filtration);
				cerr << endl;

				cerr << std::format("essential cycles (of dimension) {}", d) << endl;
				printvector<int>(essential);
				cerr << endl;
			}

			// Now, sample from the kernel, excluding unoccupied columns *or* columns
			// corresponding to cells included *after* percolating.
			int stopat = essential[this->stoppingFunction(t)];
			
			std::vector<int> leaveout;
			std::vector<int> included(stopat-start, 0);

			for (int i=stopat; i < stop; i++) leaveout.push_back(this->filtration[i]-start);
			for (int i=start; i < stopat; i++) included[i-start] = this->filtration[i]-start;

			SparseVector<RingLike> sample = arithmetic::submatrixKernelSample<RingLike>(
				this->complex->Coboundary.Matrices[d],	// complete dth coboundary matrix
				this->coefficients,							// field
				leaveout,								// rows to exclude
				this->intuniform,						// uniform random over field
				this->RNG,								// RNG
				resources,								// compute resources
				this->DEBUG					// debugging
			);

			// If we're debugging, check whether we actually sampled from the kernel.
			if (this->DEBUG) {
				// Copy and pare down the coboundary matrix.
				SparseMatrix<RingLike> cbd = this->complex->Coboundary.Matrices[d];
				for (auto i : leaveout) cbd[i].zero();
				cbd.clear_zero_row();
				cbd.compress();

				// Multiply, and check whether there's anything in the resulting std::vector;
				// there shouldn't be (i.e. it should have size 0).
				SparseVector<RingLike> outcome = arithmetic::SparseRightMultiplication<RingLike>(
					cbd,
					sample,
					this->coefficients
				);
				assert(outcome.size() < 1);
			}

			// Variable names are confusing; sorry.
			// Knock off the offset value (so the kth entry in `essential` is the index
			// of each d-cell with respect to the filtration on the d-skeleton, not the
			// entire filtration), then erase the homological percolation times that
			// are greater than the number of cells included.
			int includes = included.size();
			std::for_each(essential.begin(), essential.end(), [start](int &k) { k -= start; });
			std::erase_if(essential, [includes](int &t) { return t > includes; });

			state.cochain = sample;
			state.includes = included;
			state.essential = essential;
			state.rank = essential.size();
			state.satisfied = include.size();	// TODO satisfied of the previous cochain
			state.occupied = included.size();
			state.energy = -exclude.size();		// TODO energy of the previous cochain
			state.t = t;

			return state;
		}


		template <typename RingLike>
		InvadedCluster<RingLike>::State InvadedCluster<RingLike>::initialize(
			InvadedCluster<RingLike>::State& state
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
		InvadedCluster<RingLike>::State InvadedCluster<RingLike>::initialize(
			SparseVector<RingLike> c,
			InvadedCluster<RingLike>::State& state
		) {
			state.cochain = c;
			return state;
		}


		template <typename RingLike>
		InvadedCluster<RingLike>::InvadedCluster(
			complexes::Complex<RingLike>* complex,
			ModelParameters parameters
		) : InvadedCluster<RingLike> (
			complex,
			parameters.coefficients,
			parameters.dimension,
			parameters.stoppingFunction,
			parameters.DEBUG
		) { };

		template <typename RingLike>
		InvadedCluster<RingLike>::InvadedCluster(
			complexes::Complex<RingLike>* complex,
			Ring* R,
			int dimension,
			std::function<int(int)> stoppingFunction,
			bool DEBUG
		) : Model<RingLike,SparseVector>(R, dimension, DEBUG) {
			this->complex = complex;
			this->stoppingFunction = stoppingFunction;

			// Determine the field and build the boundary matrices for the Complex.
			this->complex->constructBoundaryMatrices(this->coefficients);

			// If the characteristic of the field is 0 (i.e. we're over the
			// rationals), then compute persistence over Z/3Z.
			if (this->coefficients->characteristic == 0) {
				Zp RR(3);
				this->coefficients = &RR;
			}

			// If the characteristic of the field is greater than 2 (i.e. we're
			// spoofing the rationals or are not in Z/2Z)...
			if (this->coefficients->characteristic > 2) {
				this->intuniform = std::uniform_int_distribution<int>(0, this->coefficients->characteristic-1);
				this->complex->constructFullBoundaryMatrix(this->coefficients);
			} else {
				this->intuniform = std::uniform_int_distribution<int>(0, 1);
				this->complex->constructFlatBoundaryMatrix();
			}
			
			std::vector<int> filtration(this->complex->size());
			std::iota(filtration.begin(), filtration.end(), 0);
			this->filtration = filtration;

			// Default set of indices to compare against.
			std::vector<int> include(this->complex->Cells[this->dimension]);
			std::iota(include.begin(), include.end(), 0);
			this->include = include;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
		};
	}
}

#endif
