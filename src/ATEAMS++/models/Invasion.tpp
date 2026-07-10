
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
		template <typename RingLike>
		ModelState<RingLike,DenseVector> Invasion<RingLike>::sample(
			int t,
			ModelState<RingLike,DenseVector>& state,
			arithmetic::ThreadOptions& options
		) {
			int d = this->dimension;

			// Shuffle indices.
			std::shuffle(this->indices.begin(), this->indices.end(), this->RNG);

			int start = this->complex->Breaks[d][0];
			int stop = this->complex->Breaks[d][1];
			int offset = (d > 0) ? this->complex->Offsets[d-1] : 0;
			
			// Fill in all the indices for cells of dimension not equal to d.
			// TODO a small optimization: we only have to do the persistence algorithm
			// through the (d+1)-dimensional cells to capture the d-dimensional events.
			// This will probably save a boatload of time!
			for (int j=0; j < this->indices.size(); j++) this->filtration[j+start] = this->indices[j]+start;

			// Persist.
			DenseVector<int> essential;

			int mod = (int)this->coefficients->characteristic;

			if (mod < 3) {
				essential = arithmetic::PHATPersistence<RingLike>(this->complex, this->filtration, this->dimension);
			} else {
				essential = arithmetic::TwistPersistence<RingLike>(this->complex, this->filtration, this->coefficients, this->dimension);
			}

			// Make sure we shift the listed indices lower.
			std::erase_if(essential, [stop, start](int x) { return !((start <= x) && (x < stop)); });

			// Now, determine at what time we hit the desired number of giant cycles
			// and return the indices of the d-cells included at or before that time.
			std::for_each(essential.begin(), essential.end(), [offset](int &k) { k -= offset; });
			std::sort(essential.begin(), essential.end()); // might be superfluous
			
			int stopAtTime = essential[this->stoppingFunction(t)];
			DenseVector<int> included(this->indices.begin(), this->indices.begin()+stopAtTime);

			state.includes = included;
			state.essential = essential;
			state.t = t;

			return state;
		}


		template <typename RingLike>
		Invasion<RingLike>::Invasion(
			complexes::Complex<RingLike>* complex,
			ModelParameters parameters
		) : Model<RingLike,DenseVector>(
			parameters.coefficients,
			parameters.dimension,
			parameters.DEBUG
		) {
			this->complex = complex;
			this->stoppingFunction = parameters.stoppingFunction;

			// Determine the field and build the boundary matrices for the Complex.
			this->complex->constructBoundaryMatrices(this->coefficients);

			int mod = (int)this->coefficients->characteristic;
			this->intuniform = std::uniform_int_distribution<int>(0,mod);

			if (mod< 3) this->complex->constructFlatBoundaryMatrix();
			else this->complex->constructFullBoundaryMatrix(this->coefficients);

			// Initialize the indices so we don't have to re-create them over and over.
			DenseVector<int> include(this->complex->Cells[this->dimension]);
			iota(begin(include), end(include), 0);
			this->indices = include;

			// Also initialize the filtration.
			DenseVector<int> filtration(this->complex->size());
			iota(begin(filtration), end(filtration), 0);
			this->filtration = filtration;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
		};

		
		template <typename RingLike>
		Invasion<RingLike>::Invasion(
			ATEAMS::complexes::Complex<RingLike>* complex,
			Ring* R,
			int dimension,
			std::function<int(int)> stoppingFunction,
			bool DEBUG
		) : Model<RingLike,DenseVector>(R, dimension, DEBUG) {
			this->complex = complex;
			this->stoppingFunction = stoppingFunction;

			// Determine the field and build the boundary matrices for the Complex.
			this->complex->constructBoundaryMatrices(this->coefficients);

			int mod = (int)this->coefficients->characteristic;
			this->intuniform = std::uniform_int_distribution<int>(0,mod);

			if (mod < 3) this->complex->constructFlatBoundaryMatrix();
			else this->complex->constructFullBoundaryMatrix(this->coefficients);

			// Initialize the indices so we don't have to re-create them over and over.
			DenseVector<int> include(this->complex->Cells[this->dimension]);
			iota(begin(include), end(include), 0);
			this->indices = include;

			// Also initialize the filtration.
			DenseVector<int> filtration(this->complex->size());
			iota(begin(filtration), end(filtration), 0);
			this->filtration = filtration;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
		};
	}
}

#endif
