
#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/InvasionPercolation.h"
#include "ATEAMS++/arithmetic/kernel.h"
#include "ATEAMS++/arithmetic/persistence.h"

#include <SparseRREF/sparse_vec.h>

#include <random>
#include <algorithm>
#include <cassert>

using namespace ATEAMS;
using namespace std;


vector<int> models::InvasionPercolation::sample(int t, arithmetic::ThreadOptions& options) {
	int d = this->parameters.dimension;
	vector<int> indices = this->indices;
	vector<int> filtration = this->filtration;

	// Shuffle indices.
	std::shuffle(indices.begin(), indices.end(), this->RNG);

	int stop = this->complex->breaks[d][0];
	int start = this->complex->breaks[d][1];
	int offset = (d > 0) ? this->complex->offsets[d-1] : 0;
	
	// Fill in all the indices for cells of dimension not equal to d.
	// TODO a small optimization: we only have to do the persistence algorithm
	// through the (d+1)-dimensional cells to capture the d-dimensional events.
	// This will probably save a boatload of time!
	for (int j=0; j < this->indices.size(); j++) filtration[j+start] = indices[j]+start;

	cout << "wh" << endl;
	printvector<int>(filtration);

	// Persist.
	vector<int> essential;

	if (this->parameters.field < 3) {
		essential = arithmetic::PHATPersistence(this->complex, filtration, this->parameters.dimension);
	} else {
		essential = arithmetic::TwistPersistence(this->complex, filtration, this->field, this->parameters.dimension);
	}

	std::erase_if(essential, [stop, start](int x) { return !((start <= x) && (x < stop)); });
	std::sort(essential.begin(), essential.end());

	// Now, determine at what time we hit the desired number of giant cycles
	// and return the indices of the d-cells included at or before that time.
	std::for_each(essential.begin(), essential.end(), [offset](int &k) { k -= offset; });
	std::sort(essential.begin(), essential.end()); // might be superfluous

	int stopAtTime = essential[this->parameters.stoppingFunction(t)];
	printvector<int>(essential);
	cout << stopAtTime << endl;

	// vector<int> included(this->indices.begin(), this->indices.end());
	// std::erase_if(included, [stopAtTime, offset](int x) { return !(x <= stopAtTime-offset); });
	// printvector<int>(included);

	vector<int> x;
	return x;
}


models::InvasionPercolation::InvasionPercolation(complexes::Complex* complex, InvasionPercolationParameters parameters)
	: field(Zp(SparseRREF::FIELD_Fp, parameters.field))
{
	this->parameters = parameters;
	this->complex = complex;

	// Determine the field and build the boundary matrices for the Complex.
	this->complex->constructBoundaryMatrices(this->field);

	if (this->parameters.field < 3) this->complex->constructFlatBoundaryMatrix();
	else this->complex->constructFullBoundaryMatrix(this->field);

	// Initialize the indices so we don't have to re-create them over and over.
	vector<int> include(this->complex->Cells[this->parameters.dimension]);
	iota(begin(include), end(include), 0);
	this->indices = include;

	// Also initialize the filtration.
	vector<int> filtration(this->complex->size());
	iota(begin(filtration), end(filtration), 0);
	this->filtration = filtration;

	// Initialize a random number generator.
	std::random_device rd;
	this->RNG = std::mt19937(rd());
	this->intuniform = std::uniform_int_distribution<int>(0,this->parameters.field);
}
