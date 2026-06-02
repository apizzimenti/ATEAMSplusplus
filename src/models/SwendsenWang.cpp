
#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/SwendsenWang.h"
#include "ATEAMS++/arithmetic/kernel.h"
#include "libraries/SparseRREF/sparse_mat.h"

#include <random>
#include <algorithm>
#include <cassert>

using namespace ATEAMS;
using namespace std;

/**
 * Samples from the conditional distribution f ~ PLGT(- | P), where P is a
 * subcomplex of the ambient complex on which f is a cocycle. 
 */
ZpVector SwendsenWang::sample(int t) {
	// Compute the temperature and probability of including particular (d-1)-cells.
	double T = this->parameters.temperatureFunction(t);
	double p = 1-exp(T);

	// If we're debugging, check that the inclusion probability is actually within
	// acceptable parameters.
	if (this->parameters._DEBUG) {
		cerr << std::format("verifying 0 ≤ {} ≤ 1", p) << endl;
		assert((0 <= p) && (p <= 1));
	}

	// Perform independent Bernoulli trials over each (d-1)-cell to decide whether
	// it gets included in the computation. We do this by performing a (sparse!) matrix
	// multiplication, then checking which of the entries are zero.
	size_t d = this->parameters.dimension;

	ZpVector coefficients = sparse_mat_dot_sparse_vec<data_t, index_t>(
		this->complex->Coboundary.Matrices[d],
		this->state.cochain,
		this->field
	);
	coefficients.compress();

	// The zero entries are those that don't show up in the indices, so we
	// automatically exclude these.
	int N = this->complex->Cells[d];
	set<size_t> exclude;
	
	if (coefficients.size() > 0) {
		for (int e=0; e < coefficients.size(); e++) exclude.insert((size_t)e);
	}

	// We then iterate over the remaining indices and flip a weighted coin (i.e.
	// perform a density-p Bernoulli trial) to decide whether to keep that
	// particular row (corresponding to a d-cell on which the current cochain
	// evaluates to 0).
	for (int i=0; i<N; i++) {
		if (this->unituniform(this->RNG) > p) exclude.insert(i);
	}

	// Now, sample from the kernel.
	ZpVector sample = submatrixKernelSample(
		this->complex->Coboundary.Matrices[d],	// full coboundary matrix
		this->field,							// field
		exclude,								// rows to exclude
		this->intuniform,						// uniform random over field
		this->RNG,								// RNG
		this->parameters._DEBUG					// debugging
	);

	// If we're debugging, check whether we actually sampled from the kernel.
	if (this->parameters._DEBUG) {
		// Copy and pare down the coboundary matrix.
		ZpMatrix cbd = this->complex->Coboundary.Matrices[d];
		for (auto i : exclude) cbd[i].zero();
		cbd.clear_zero_row();
		cbd.compress();

		// Multiply, and check whether there's anything in the resulting vector;
		// there shouldn't be (i.e. it should have size 0).
		ZpVector outcome = sparse_mat_dot_sparse_vec<data_t, index_t>(
			cbd,
			sample,
			this->field
		);

		assert(outcome.size() < 1);
	}

	this->state.cochain = sample;
	
	return sample;
}

/**
 * Creates an initial cochain of randoms (in the appropriate range).
 */
void SwendsenWang::initial() {
	size_t dimension = this->parameters.dimension-1;
	int N = this->complex->Cells[dimension];

	ZpVector cochain;
	for (int i=0; i < N; i++) cochain.push_back(
		(index_t)i, (data_t)this->intuniform(this->RNG)
	);

	cochain.compress();
	this->state.cochain = cochain;
}

/**
 * Sets an initial cochain from the provided vector.
 */
void SwendsenWang::initial(ZpVector c) {
	this->state.cochain = c;
}

/**
 * Default constructor.
 */
SwendsenWang::SwendsenWang(Complex* complex, SwendsenWangParameters parameters)
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
}
