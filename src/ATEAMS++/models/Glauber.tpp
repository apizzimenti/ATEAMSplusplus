
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

namespace ATEAMS {
	namespace models {
		template <typename T>
		models::ModelState<T,SparseVector> models::Glauber<T>::sample(
			int t,
			models::ModelState<T,SparseVector>& state,
			arithmetic::ThreadOptions& options
		) {
			// Compute the temperature and probability of including particular (d-1)-cells.
			double temp = this->parameters.temperatureFunction(t);
			
			// Sample a (d-1)-cell at uniform random and change its spin.
			index_t flippedIndex = this->indexuniform(this->RNG);
			SparseVector<T> flipped(state.cochain);

			// Overwrite the entry? Not sure if this actually overwrites it or not.
			flipped.push_back(
				(index_t)flippedIndex,
				(T)this->intuniform(this->RNG)
			);

			// Compute the energies of the cochains.
			// TODO this is gonna be really slow, so we should work out a better way to
			// compute these.
			int original = statistics::energy(this->complex, state.cochain, this->field, this->parameters.dimension);
			int updated = statistics::energy(this->complex, flipped, this->field, this->parameters.dimension);

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

			return state;
		}


		template <typename T>
		models::ModelState<T,SparseVector> models::Glauber<T>::initialize(models::ModelState<T,SparseVector>& state) {
			size_t dimension = this->parameters.dimension-1;
			int N = this->complex->Cells[dimension];

			SparseVector<T> cochain;
			for (int i=0; i < N; i++) cochain.push_back(
				(index_t)i, (T)this->intuniform(this->RNG)
			);

			cochain.compress();
			state.cochain = cochain;

			return state;
		}


		template <typename T>
		models::ModelState<T,SparseVector> models::Glauber<T>::initialize(
			SparseVector<T> c,
			models::ModelState<T,SparseVector>& state
		) {
			state.cochain = c;
			return state;
		}

		template <typename T>
		models::Glauber<T>::Glauber(
			complexes::Complex<T>* complex,
			models::ModelParameters parameters
		)
			: field(parameters.field > 0 ? Field(SparseRREF::FIELD_Fp, parameters.field) : Field(SparseRREF::FIELD_QQ))
		{
			this->parameters = parameters;
			this->complex = complex;

			// Determine the field and build the boundary matrices for the Complex.
			this->complex->constructBoundaryMatrices(this->field);

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->unituniform = std::uniform_real_distribution<double>(0,1);
			this->intuniform = std::uniform_int_distribution<int>(0,this->parameters.field > 0 ? this->parameters.field : 1);
			this->indexuniform = std::uniform_int_distribution<int>(0, this->complex->Cells[this->parameters.dimension-1]);
		}
	}
}

#endif
