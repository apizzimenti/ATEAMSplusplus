
#include "ATEAMS++/models/Bernoulli.h"
#include "ATEAMS++/arithmetic/persistence.h"

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
	this->state.included = include;

	// Construct the filtration.
	vector<int> filtration(this->complex->offsets[d+1], 0);
	int stop = this->complex->breaks[d][0];
	int start = this->complex->breaks[d][1];
	int offset = (d > 0) ? this->complex->offsets[d-1] : 0;
	
	// Fill in all the indices for cells of dimension not equal to d.
	// TODO a small optimization: we only have to do the persistence algorithm
	// through the (d+1)-dimensional cells to capture the d-dimensional events.
	// This will probably save a boatload of time!
	for (int j=0; j < stop; j++) filtration[j] = j;
	for (int j=start; j < this->complex->offsets[d+1]; j++) filtration[j] = j;

	for (int j=0; j < include.size(); j++) filtration[j+offset] = include[j]+offset;
	for (int j=0; j < exclude.size(); j++) filtration[j+offset+include.size()] = exclude[j]+offset;

	// TODO finish this debugging thing
	// if (this->parameters.DEBUG && this->complex->periodic) {
	// 	std::cerr << "verifying whether each d-cell is shared by the correct number of (d+1)-cells..." << std::endl;
	// }

	// Now, compute the persistence times, then filter over them to capture only
	// the ones within the right window.
	vector<int> essential = arithmetic::PHATPersistence(this->complex, filtration, d);
	std::erase_if(essential, [stop, included](int t) { return !((stop <= t) && (t < stop+included)); });
	std::sort(essential.begin(), essential.end());

	if (this->parameters.DEBUG) {
		std::cerr << "outputting filtration to stdout..." << std::endl;
		printvector(filtration);

		std::cerr << "outputting essential cycles to stdout..." << std::endl;
		printvector(essential);
	}

	// Knock off the offset value, set state, and done.
	std::for_each(essential.begin(), essential.end(), [offset](int &k) { k -= offset; });
	this->state.essential = essential;
	this->state.rank = essential.size();

	return include;
}


models::Bernoulli::Bernoulli(complexes::Complex* complex, BernoulliParameters parameters) 
	: field(Zp(SparseRREF::FIELD_Fp, 2)) 
{
	this->parameters = parameters;
	this->complex = complex;

	// From those boundary matrices, construct the "full" (i.e. up to dimension d+1)
	// PHAT boundary matrix.
	this->complex->constructFlatBoundaryMatrix();
	this->complex->constructFullBoundaryMatrix(this->field);


	// Initialize a random number generator.
	std::random_device rd;
	this->RNG = std::mt19937(rd());
	this->unituniform = std::uniform_real_distribution<double>(0,1);
}
