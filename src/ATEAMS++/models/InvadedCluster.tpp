
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
			arithmetic::ComputeOptions& options
		) {
			int d = this->dimension;

			// Determine which (d-1)-cells are "satisfiable" (i.e. those on which the
			// current (d-1)-cochain vanishes).
			SparseVector<RingLike> satisfiable = arithmetic::SparseRightMultiplication<RingLike>(
				this->complex->Coboundary.Matrices[d],
				state.cochain,
				this->coefficients
			);

			// The zero entries are those that don't show up in the indices, so we
			// automatically exclude these.
			int N = this->complex->Cells[d];
			std::set<size_t> exclude;
			std::vector<size_t> include;
			
			if (satisfiable.size() > 0) {
				for (int e=0; e < satisfiable.size(); e++) exclude.insert((size_t)e);
			}

			for (int i=0; i < N; i++) {
				if (!exclude.contains(i)) include.push_back(i);
			}

			std::shuffle(include.begin(), include.end(), this->RNG);

			// Construct the filtration.
			std::vector<size_t> inside(include.begin(), include.end());
			inside.insert(inside.end(), exclude.begin(), exclude.end());

			int stop = this->complex->Breaks[d][1];
			int start = this->complex->Breaks[d][0];
			
			// Fill in all the indices for cells of dimension not equal to d.
			// TODO a small optimization: we only have to do the persistence algorithm
			// through the (d+1)-dimensional cells to capture the d-dimensional events.
			// This will probably save a boatload of time!
			for (int j=0; j < inside.size(); j++) this->filtration[j+start] = inside[j]+start;

			// Persist.
			std::vector<int> essential = topology::persistence<RingLike>(
				this->complex,
				this->filtration,
				this->coefficients,
				this->dimension,
				options
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
			
			std::set<size_t> leaveout;
			std::vector<int> included(stopat-start, 0);

			for (int i=stopat; i < stop; i++) leaveout.insert(this->filtration[i]-start);
			for (int i=start; i < stopat; i++) included[i-start] = this->filtration[i]-start;

			SparseVector<RingLike> sample = arithmetic::submatrixKernelSample<RingLike>(
				this->complex->Coboundary.Matrices[d],	// complete dth coboundary matrix
				this->coefficients,							// field
				leaveout,								// rows to exclude
				this->intuniform,						// uniform random over field
				this->RNG,								// RNG
				options,								// compute options
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

			state.cochain = sample;
			state.includes = included;
			state.essential = essential;
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
		) : Model<RingLike,SparseVector>(
			parameters.coefficients,
			parameters.dimension,
			parameters.DEBUG
		) {
			this->complex = complex;

			// Determine the field and build the boundary matrices for the Complex.
			this->complex->constructBoundaryMatrices(this->coefficients);
			this->stoppingFunction = parameters.stoppingFunction;

			// Set up the arithmetic depending on the coefficient field. Currently,
			// we use Z/2Z arithmetic when computing over the rationals.
			if (this->coefficients->characteristic > 2) {
				this->intuniform = std::uniform_int_distribution<int>(0, this->coefficients->characteristic);
				this->complex->constructFullBoundaryMatrix(this->coefficients);
			} else {
				this->intuniform = std::uniform_int_distribution<int>(0, 2);
				this->complex->constructFlatBoundaryMatrix();
			}

			std::vector<int> filtration(this->complex->size());
			std::iota(filtration.begin(), filtration.end(), 0);
			this->filtration = filtration;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
		};

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

			// Set up the arithmetic depending on the coefficient field. Currently,
			// we use Z/2Z arithmetic when computing over the rationals.
			if (this->coefficients->characteristic > 2) {
				this->intuniform = std::uniform_int_distribution<int>(0, this->coefficients->characteristic);
				this->complex->constructFullBoundaryMatrix(this->coefficients);
			} else {
				this->intuniform = std::uniform_int_distribution<int>(0, 2);
				this->complex->constructFlatBoundaryMatrix();
			}

			std::vector<int> filtration(this->complex->size());
			std::iota(filtration.begin(), filtration.end(), 0);
			this->filtration = filtration;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
		};
	}
}

#endif
