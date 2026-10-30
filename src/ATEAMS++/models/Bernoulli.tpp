
#ifndef ATEAMS_MODELS_BERNOULLI_T
#define ATEAMS_MODELS_BERNOULLI_T

#ifndef ATEAMS_MODELS_BERNOULLI_H
#error __FILE__ should only be included from models/Bernoulli.h.
#endif

#include "ATEAMS++/models/Bernoulli.h"
#include "ATEAMS++/topology/persistence.h"

#include <random>
#include <cassert>

namespace ATEAMS::models {
	Bernoulli::State Bernoulli::sample(
		int t,
		Bernoulli::State& state,
		arithmetic::ComputeResources<RingType>& resources
	) {
		// If we're debugging, check that the inclusion probability is a probability
		// to begin with.
		double p = this->p;
		int d = this->dimension;
		
		if (this->DEBUG) {
			cerr << std::format("verifying 0 ≤ {} ≤ 1", p) << endl;
			assert((0 <= p) && (p <= 1));
		}

		// Flip a density-p coin for each cell.
		vector<double>& probabilities = this->probabilities;

		resources.random->fill<std::mt19937,std::uniform_real_distribution<double>>(
			this->RNG, this->unituniform, probabilities
		);

		// Check whether the coinflips are less than p, organizing in/excluded
		// indices. Requires one pass over all cells.
		int included = 0,
			excluded = 0;

		vector<int>& include = this->include;
		vector<int>& exclude = this->exclude;

		for (size_t t=0; t < this->complex->Cells[d]; t++) {
			if (probabilities[t] < p) {
				include[included] = t;
				included++;
			} else {
				exclude[excluded] = t;
				excluded++;
			}
		}

		// Construct the filtration.
		int stop = this->complex->Breaks[d][0],
			start = this->complex->Breaks[d][1],
			offset = (d > 0) ? this->complex->Offsets[d-1] : 0;
		
		// Fill in all the indices for cells of dimension not equal to d.
		for (int j=0; j < stop; j++) this->filtration[j] = j;
		for (int j=start; j < this->complex->Offsets[d+1]; j++) this->filtration[j] = j;
		
		// Shuffle the included indices.
		std::shuffle(include.begin(), include.begin()+included, this->RNG);

		// Build out the middle part of the filtration.
		for (int j=0; j < included; j++) this->filtration[j+offset] = include[j]+offset;
		for (int j=0; j < excluded; j++) this->filtration[j+offset+included] = exclude[j]+offset;
		
		// Now, compute the persistence times, then filter over them to capture only
		// the ones within the right window.
		vector<int> essential = topology::persistence::persistence<RingType>(this->complex, this->filtration, this->coefficients, d, resources);

		if (this->DEBUG) {
			std::cerr << "outputting filtration to stdout..." << std::endl;
			printvector<int>(this->filtration);

			std::cerr << "outputting essential cycles to stdout..." << std::endl;
			printvector<int>(essential);
		}

		// Knock off the offset value (so the kth entry in `essential` is the index
		// of each d-cell with respect to the filtration on the d-skeleton, not the
		// entire filtration), then erase the homological percolation times that
		// are greater than the number of cells included.
		std::for_each(essential.begin(), essential.end(), [offset](int &k) { k -= offset; });
		std::erase_if(essential, [included](int &t) { return t > included; });

		state.includes = include; 
		state.essential = essential;
		state.rank = essential.size();
		state.t = t;

		return state;
	}

	Bernoulli::Bernoulli(
		complexes::Complex<RingType>* complex,
		ModelParameters parameters
	) : Bernoulli(
		complex,
		parameters.dimension,
		parameters.p,
		parameters.DEBUG
	) { };

	Bernoulli::Bernoulli(
		complexes::Complex<RingType>* complex,
		int dimension,
		double p,
		bool DEBUG
	) : Model<RingType,DenseVector>(new RingType, dimension, DEBUG) {
		this->complex = complex;
		this->p = p;

		// From those boundary matrices, construct the "full" (i.e. up to dimension d+1)
		// PHAT boundary matrix.
		this->complex->constructFlatBoundaryMatrix();

		// Resuable containers for filtrations, probabilities, included/excluded vectors.
		this->filtration = vector<int>(this->complex->Offsets[this->dimension+1]);
		this->probabilities = vector<double>(this->complex->Cells[dimension]);
		this->include = vector<int>(this->complex->Cells[dimension]);
		this->exclude = vector<int>(this->complex->Cells[dimension]);

		// Initialize a random number generator.
		std::random_device rd;
		this->RNG = std::mt19937(rd());
		this->unituniform = std::uniform_real_distribution<double>(0,1);
	};
}

#endif
