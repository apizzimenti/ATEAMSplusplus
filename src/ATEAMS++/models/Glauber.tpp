
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

#include <random>
#include <algorithm>
#include <cassert>

namespace ATEAMS::models {
		
	template <typename RingLike>
	Glauber<RingLike>::State Glauber<RingLike>::sample(
		int t,
		Glauber<RingLike>::State& state,
		arithmetic::ComputeOptions<RingLike>& options
	) {
		// Compute the temperature and probability of including particular (d-1)-cells.
		double temp = this->temperatureFunction(t);
		
		// Sample a (d-1)-cell at uniform random and change its spin.
		INDEX flippedIndex = this->indexuniform(this->RNG);
		SparseVector<RingLike> flipped(state.cochain);

		// Overwrite the entry? Not sure if this actually overwrites it or not.
		flipped.push_back(
			(INDEX)flippedIndex,
			(typename RingLike::dtype)this->intuniform(this->RNG)
		);

		// Compute the energies of the cochains.
		// TODO this is gonna be really slow, so we should work out a better way to
		// compute these.
		int original = statistics::energy(this->complex, state.cochain, this->coefficients, this->dimension);
		int updated = statistics::energy(this->complex, flipped, this->coefficients, this->dimension);

		// Comptue the probability of accepting the new cochain; perform the Metropolis
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

		return state;
	}


	template <typename RingLike>
	Glauber<RingLike>::State Glauber<RingLike>::initialize(
		SparseVector<RingLike> c,
		Glauber<RingLike>::State& state
	) {
		state.cochain = c;
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
		this->indexuniform = std::uniform_int_distribution<int>(0, this->complex->Cells[this->dimension-1]);

		this->intuniform = std::uniform_int_distribution<int>(
			0,
			this->coefficients->characteristic > 0 ? this->coefficients->characteristic : 2
		);
	};

	template <typename RingLike>
	Glauber<RingLike>::Glauber(
		complexes::Complex<RingLike>* complex,
		Ring* R,
		int dimension,
		bool DEBUG
	) : Model<RingLike,SparseVector>(R, dimension, DEBUG) {
		this->temperatureFunction = statistics::selfdual(this->coefficients);
		this->complex = complex;

		// Determine the field and build the boundary matrices for the Complex.
		this->complex->constructBoundaryMatrices(this->coefficients);

		// Initialize a random number generator.
		std::random_device rd;
		this->RNG = std::mt19937(rd());
		this->unituniform = std::uniform_real_distribution<double>(0,1);
		this->indexuniform = std::uniform_int_distribution<int>(0, this->complex->Cells[this->dimension-1]);

		this->intuniform = std::uniform_int_distribution<int>(
			0,
			this->coefficients->characteristic > 0 ? this->coefficients->characteristic : 2
		);
	};


	template <typename RingLike>
	Glauber<RingLike>::Glauber(
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

		// Initialize a random number generator.
		std::random_device rd;
		this->RNG = std::mt19937(rd());
		this->unituniform = std::uniform_real_distribution<double>(0,1);
		this->indexuniform = std::uniform_int_distribution<int>(0, this->complex->Cells[this->dimension-1]);

		this->intuniform = std::uniform_int_distribution<int>(
			0,
			this->coefficients->characteristic > 0 ? this->coefficients->characteristic : 2
		);
	};
}

#endif
