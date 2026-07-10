
#ifndef ATEAMS_MODELS_BERNOULLI_T
#define ATEAMS_MODELS_BERNOULLI_T

#ifndef ATEAMS_MODELS_BERNOULLI_H
#error __FILE__ should only be included from models/Bernoulli.h.
#endif

#include "ATEAMS++/models/Bernoulli.h"
#include "ATEAMS++/topology/persistence.h"

#include <random>
#include <cassert>

namespace ATEAMS {
	namespace models {
		ModelState<Z2,DenseVector> Bernoulli::sample(
			int t,
			ModelState<Z2,DenseVector>& state,
			arithmetic::ThreadOptions& options
		) {
			// If we're debugging, check that the inclusion probability is a probability
			// to begin with.
			double p = this->p;
			int d = this->dimension;
			
			if (this->DEBUG) {
				cerr << std::format("verifying 0 ≤ {} ≤ 1", p) << endl;
				assert((0 <= p) && (p <= 1));
			}

			// Then decide which cells we're including/excluding.
			vector<int> include;
			vector<int> exclude;

			for (int i=0; i < this->complex->Cells[d]; i++) {
				if (this->unituniform(this->RNG) < p) {
					include.push_back(i);
				} else {
					exclude.push_back(i);
				}
			}

			int included = include.size();
			state.includes = include;

			// Construct the filtration.
			int stop = this->complex->Breaks[d][0];
			int start = this->complex->Breaks[d][1];
			int offset = (d > 0) ? this->complex->Offsets[d-1] : 0;
			
			// Fill in all the indices for cells of dimension not equal to d.
			for (int j=0; j < stop; j++) this->filtration[j] = j;
			for (int j=start; j < this->complex->Offsets[d+1]; j++) this->filtration[j] = j;

			// Shuffle the included indices.
			std::shuffle(include.begin(), include.end(), this->RNG);
			
			for (int j=0; j < include.size(); j++) this->filtration[j+offset] = include[j]+offset;
			for (int j=0; j < exclude.size(); j++) this->filtration[j+offset+include.size()] = exclude[j]+offset;

			// Now, compute the persistence times, then filter over them to capture only
			// the ones within the right window.
			vector<int> essential = topology::persistence<Z2>(this->complex, this->filtration, this->coefficients, d);
			std::erase_if(essential, [stop, included](int t) { return !((stop <= t) && (t < stop+included)); });

			if (this->DEBUG) {
				std::cerr << "outputting filtration to stdout..." << std::endl;
				printvector<int>(this->filtration);

				std::cerr << "outputting essential cycles to stdout..." << std::endl;
				printvector<int>(essential);
			}

			// Knock off the offset value, set state, and done.
			std::for_each(essential.begin(), essential.end(), [offset](int &k) { k -= offset; });
			state.essential = essential;
			state.rank = essential.size();
			state.t = t;

			return state;
		}

		Bernoulli::Bernoulli(
			complexes::Complex<Z2>* complex,
			ModelParameters parameters
		) : Model<Z2,DenseVector>(
			new Z2,
			parameters.dimension,
			parameters.DEBUG
		) {
			this->complex = complex;
			this->p = parameters.p;

			// From those boundary matrices, construct the "full" (i.e. up to dimension d+1)
			// PHAT boundary matrix.
			this->complex->constructFlatBoundaryMatrix();
			this->complex->constructFullBoundaryMatrix(this->coefficients);

			// Default filtration.
			vector<int> filtration(this->complex->Offsets[this->dimension+1], 0);
			this->filtration = filtration;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->unituniform = std::uniform_real_distribution<double>(0,1);
		};

		Bernoulli::Bernoulli(
			complexes::Complex<Z2>* complex,
			int dimension,
			double p,
			bool DEBUG
		) : Model<Z2,DenseVector>(new Z2, dimension, DEBUG) {
			this->complex = complex;

			// From those boundary matrices, construct the "full" (i.e. up to dimension d+1)
			// PHAT boundary matrix.
			this->complex->constructFlatBoundaryMatrix();
			this->complex->constructFullBoundaryMatrix(this->coefficients);

			// Default filtration.
			vector<int> filtration(this->complex->Offsets[this->dimension+1], 0);
			this->filtration = filtration;

			// Initialize a random number generator.
			std::random_device rd;
			this->RNG = std::mt19937(rd());
			this->unituniform = std::uniform_real_distribution<double>(0,1);
		};
	}
}

#endif
