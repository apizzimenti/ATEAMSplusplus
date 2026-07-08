
#ifndef ATEAMS_MODELS_INVASION_T
#define ATEAMS_MODELS_INVASION_T

#ifndef ATEAMS_MODELS_INVASION_H
#error __FILE__ should only be included from models/Invasion.h.
#endif

#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Invasion.h"
#include "ATEAMS++/arithmetic/kernel.h"
#include "ATEAMS++/arithmetic/persistence.h"

#include <SparseRREF/sparse_vec.h>

#include <random>
#include <algorithm>
#include <cassert>

using namespace std;

namespace ATEAMS {
	namespace models {
		ModelState<ff,DenseVector> Invasion::sample(
			int t,
			ModelState<ff,DenseVector>& state,
			arithmetic::ThreadOptions& options
		) {
			int d = this->parameters.dimension;

			// Shuffle indices.
			std::shuffle(this->indices.begin(), this->indices.end(), this->RNG);

			int start = this->complex->breaks[d][0];
			int stop = this->complex->breaks[d][1];
			int offset = (d > 0) ? this->complex->offsets[d-1] : 0;
			
			// Fill in all the indices for cells of dimension not equal to d.
			// TODO a small optimization: we only have to do the persistence algorithm
			// through the (d+1)-dimensional cells to capture the d-dimensional events.
			// This will probably save a boatload of time!
			for (int j=0; j < this->indices.size(); j++) this->filtration[j+start] = this->indices[j]+start;

			// Persist.
			DenseVector<int> essential;

			if (this->parameters.field < 3) {
				essential = arithmetic::PHATPersistence(this->complex, this->filtration, this->parameters.dimension);
			} else {
				essential = arithmetic::TwistPersistence(this->complex, this->filtration, this->field, this->parameters.dimension);
			}

			// Make sure we shift the listed indices lower.
			std::erase_if(essential, [stop, start](int x) { return !((start <= x) && (x < stop)); });

			// Now, determine at what time we hit the desired number of giant cycles
			// and return the indices of the d-cells included at or before that time.
			std::for_each(essential.begin(), essential.end(), [offset](int &k) { k -= offset; });
			std::sort(essential.begin(), essential.end()); // might be superfluous
			
			int stopAtTime = essential[this->parameters.stoppingFunction(t)];
			DenseVector<int> included(this->indices.begin(), this->indices.begin()+stopAtTime);

			state.includes = included;
			state.essential = essential;
			state.t = t;

			return state;
		}


		Invasion::Invasion(
			complexes::Complex<ff>* complex,
			ModelParameters parameters
		) : Model<ff,DenseVector>( // parent constructor; initializes the field.
			parameters.field > 0 ?	
				Field(SparseRREF::FIELD_Fp, parameters.field) :
				Field(SparseRREF::FIELD_QQ),
			"Invasion"
		) {
			this->parameters = parameters;
			this->complex = complex;

			// Determine the field and build the boundary matrices for the Complex.
			this->complex->constructBoundaryMatrices(this->field);

			if (this->parameters.field < 3) this->complex->constructFlatBoundaryMatrix();
			else this->complex->constructFullBoundaryMatrix(this->field);

			// Initialize the indices so we don't have to re-create them over and over.
			DenseVector<int> include(this->complex->Cells[this->parameters.dimension]);
			iota(begin(include), end(include), 0);
			this->indices = include;

			// Also initialize the filtration.
			DenseVector<int> filtration(this->complex->size());
			iota(begin(filtration), end(filtration), 0);
			this->filtration = filtration;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->intuniform = std::uniform_int_distribution<int>(0,this->parameters.field);
		}
	}
}

#endif
