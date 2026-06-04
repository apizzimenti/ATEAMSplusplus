
#include "ATEAMS++/models/Bernoulli.h"

#include <random>
#include <cassert>

using namespace ATEAMS;
using namespace std;

vector<int> models::Bernoulli::sample(int t, arithmetic::ThreadOptions& options) {
	// If we're debugging, check that the inclusion probability is a probability
	// to begin with.
	double p = this->parameters.p;
	int d = this->parameters.dimension;
	
	if (this->parameters.DEBUG) {
		cerr << std::format("verifying 0 ≤ {} ≤ 1", p) << endl;
		assert((0 <= p) && (p <= 1));
	}

	// Then decide which cells we're including.
	vector<int> s;

	for (int i=0; i < this->complex->Cells[d]; i++) {
		if (this->unituniform(this->RNG) < p) s.push_back(i);
	}

	this->state.included = s;

	// Now, compute the persistent homology.

	return s;
}


models::Bernoulli::Bernoulli(complexes::Complex* complex, BernoulliParameters parameters) 
	: field(Zp(SparseRREF::FIELD_Fp, 2)) 
{
	this->parameters = parameters;
	this->complex = complex;

	// Determine the field and build the boundary matrices for the Complex.
	this->complex->constructBoundaryMatrices(this->field);

	// From those boundary matrices, construct the "full" (i.e. up to dimension d+1)
	// PHAT boundary matrix.
	

	// Initialize a random number generator.
	std::random_device rd;
	this->RNG = std::mt19937(rd());
	this->unituniform = std::uniform_real_distribution<double>(0,1);
}
