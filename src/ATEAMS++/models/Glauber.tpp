
#ifndef ATEAMS_MODELS_GLAUBER_T
#define ATEAMS_MODELS_GLAUBER_T

#ifndef ATEAMS_MODELS_GLAUBER_H
#error __FILE__ should only be included from models/Glauber.h.
#endif

#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Glauber.h"
#include "ATEAMS++/statistics/observables.h"

#include <SparseRREF/sparse_type.h>
#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_vec.h>

#include <random>
#include <algorithm>
#include <cassert>

namespace ATEAMS::models {

	template <typename RingLike>
	INDEX entryExistsAtIndex(SparseVector<RingLike> &u, INDEX index) {
		int L = 0, R = u.size()-1, m;

		while (L <= R) {
			m = L + (R-L)/2;

			if (u.indices[m] == index) return (INDEX)m;
			else if (u.indices[m] < index) L = m+1;
			else R = m-1;
		}

		return (INDEX)(u.size() + 1);
	}
		
	template <typename RingLike>
	Glauber<RingLike>::State Glauber<RingLike>::sample(
		int t,
		Glauber<RingLike>::State& state,
		arithmetic::ComputeResources<RingLike>& resources
	) {
		// Compute the temperature and probability of including particular (d-1)-cells.
		double temp = this->temperatureFunction(t);
		
		// Sample a (d-1)-cell at uniform random and change its spin.
		INDEX flippedIndex = this->indexuniform(this->RNG);
		SparseVector<RingLike> flipped(state.cochain);

		INDEX exists = entryExistsAtIndex<RingLike>(flipped, flippedIndex);

		if (exists < flipped.size()) {
			flipped[exists] = (typename RingLike::dtype)this->intuniform(this->RNG);
		} else {
			flipped.push_back(
				(INDEX)flippedIndex,
				(typename RingLike::dtype)this->intuniform(this->RNG)
			);
		}

		flipped.compress();

		// From the boundary matrix, check the indices of the d-cells that have
		// the chosen (d-1)-cell as a face; we only need to compute the difference
		// in the energy on those cells.
		SparseVector<RingLike> coboundary = this->complex->Boundary.Matrices[this->dimension][flippedIndex];
		int originalNonzero = 0, flippedNonzero = 0;

		// #pragma omp parallel for
		for (int i=0; i < coboundary.size(); i++) {
			int index = coboundary(i);
			SparseVector<RingLike> boundary = this->complex->Coboundary.Matrices[this->dimension][index];

			originalNonzero += (int)(arithmetic::SparseVectorDot<RingLike>(state.cochain, boundary, this->coefficients) > 0);
			flippedNonzero += (int)(arithmetic::SparseVectorDot<RingLike>(flipped, boundary, this->coefficients) > 0);
		}

		int original = state.energy;
		int updated = state.energy + (originalNonzero-flippedNonzero);

		// Compute the probability of accepting the new cochain; perform the Metropolis
		// step.
		double random = this->unituniform(this->RNG);
		double diff = min(1.0, exp(temp*(original-updated)));

		if (random < diff) {
			state.cochain = flipped;
			state.energy = updated;
		} else {
			state.energy = original;
		}

		state.t = t;
		
		return state;
	}


	template <typename RingLike>
	Glauber<RingLike>::State Glauber<RingLike>::initialize(Glauber<RingLike>::State& state) {
		size_t dimension = this->dimension-1;
		int N = this->complex->Cells[dimension];

		SparseVector<RingLike> cochain;
		for (int i=0; i < N; i++) cochain.push_back(
			(INDEX)i, (typename RingLike::dtype)this->intuniform(this->RNG)
		);

		cochain.compress();
		state.cochain = cochain;
		state.energy = statistics::energy(this->complex, state.cochain, this->coefficients, this->dimension);

		return state;
	}


	template <typename RingLike>
	Glauber<RingLike>::State Glauber<RingLike>::initialize(
		SparseVector<RingLike> c,
		Glauber<RingLike>::State& state
	) {
		state.cochain = c;
		state.energy = statistics::energy(this->complex, state.cochain, this->coefficients, this->dimension);
		return state;
	};

	template <typename RingLike>
	Glauber<RingLike>::Glauber(
		complexes::Complex<RingLike>* complex,
		Ring* R,
		int dimension,
		std::function<double(int)> temperatureFunction,
		bool DEBUG
	) : Model<RingLike,SparseVector>(R, dimension, DEBUG) {
		this->temperatureFunction = temperatureFunction;
		this->complex = complex;

		// Determine the field and build the boundary matrices for the Complex.
		this->complex->constructBoundaryMatrices(this->coefficients);

		// Initialize a random number generator.
		std::random_device rd;
		this->RNG = std::mt19937(rd());
		this->unituniform = std::uniform_real_distribution<double>(0,1);
		this->indexuniform = std::uniform_int_distribution<int>(0, this->complex->Cells[this->dimension-1]-1);

		this->intuniform = std::uniform_int_distribution<int>(
			0,
			this->coefficients->characteristic > 0 ? this->coefficients->characteristic-1 : 1
		);
	};

	template <typename RingLike>
	Glauber<RingLike>::Glauber(
		complexes::Complex<RingLike>* complex,
		Ring* R,
		int dimension,
		bool DEBUG
	) : Glauber<RingLike>(
		complex,
		R,
		dimension,
		statistics::selfdual,
		DEBUG
	) { };


	template <typename RingLike>
	Glauber<RingLike>::Glauber(
		complexes::Complex<RingLike>* complex,
		ModelParameters parameters
	) : Glauber<RingLike>(
		complex,
		parameters.coefficients,
		parameters.dimension,
		parameters.temperatureFunction,
		parameters.DEBUG
	) { };
}

#endif
