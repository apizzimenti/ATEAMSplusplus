
#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/Glauber.h"
#include "ATEAMS++/statistics/observables.h"
#include <SparseRREF/sparse_type.h>
#include <SparseRREF/sparse_mat.h>

#include <random>
#include <algorithm>
#include <cassert>

using namespace ATEAMS;
using namespace std;


ZpVector models::Glauber::sample(int t, arithmetic::ThreadOptions& options) {
	// Compute the temperature and probability of including particular (d-1)-cells.
	double T = this->parameters.temperatureFunction(t);
	
	// Sample a (d-1)-cell at uniform random and change its spin.
	index_t flippedIndex = this->indexuniform(this->RNG);
	ZpVector flipped(this->state.cochain);

	// Overwrite the entry? Not sure if this actually overwrites it or not.
	flipped.push_back(
		(index_t)flippedIndex,
		(data_t)this->intuniform(this->RNG)
	);

	// Compute the energies of the cochains.
	// TODO this is gonna be really slow, so we should work out a better way to
	// compute these.
	int original = statistics::energy(this->complex, this->state.cochain, this->field, this->parameters.dimension);
	int updated = statistics::energy(this->complex, flipped, this->field, this->parameters.dimension);

	// Comptue the probability of accepting the new cochain; perform the Metropolis
	// step.
	double random = this->unituniform(this->RNG);
	double diff = min(1.0, exp(T*(original-updated)));

	if (random < diff) {
		this->state.cochain = flipped;
		this->state.energy = updated;
	} else {
		this->state.energy = original;
	}

	return this->state.cochain;
}


void models::Glauber::initialize() {
	size_t dimension = this->parameters.dimension-1;
	int N = this->complex->Cells[dimension];

	ZpVector cochain;
	for (int i=0; i < N; i++) cochain.push_back(
		(index_t)i, (data_t)this->intuniform(this->RNG)
	);

	cochain.compress();
	this->state.cochain = cochain;
}


void models::Glauber::initialize(ZpVector c) {
	this->state.cochain = c;
}

models::Glauber::Glauber(complexes::Complex* complex, GlauberParameters parameters)
	: field(Zp(SparseRREF::FIELD_Fp, parameters.field))
{
	this->parameters = parameters;
	this->complex = complex;

	// Determine the field and build the boundary matrices for the Complex.
	this->complex->constructBoundaryMatrices(this->field);

	// Initialize a random number generator.
	std::random_device rd;
	this->RNG = std::mt19937(rd());
	this->unituniform = std::uniform_real_distribution<double>(0,1);
	this->intuniform = std::uniform_int_distribution<int>(0,this->parameters.field);
	this->indexuniform = std::uniform_int_distribution<int>(0, this->complex->Cells[this->parameters.dimension-1]);
}
