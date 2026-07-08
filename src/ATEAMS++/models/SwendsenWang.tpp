
#ifndef ATEAMS_MODELS_SWENDSENWANG_T
#define ATEAMS_MODELS_SWENDSENWANG_T

#ifndef ATEAMS_MODELS_SWENDSENWANG_H
#error __FILE__ should only be included from models/SwendsenWang.h.
#endif

#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/SwendsenWang.h"
#include "ATEAMS++/arithmetic/kernel.h"

#include <SparseRREF/sparse_mat.h>

#include <random>
#include <algorithm>
#include <cassert>

using namespace ATEAMS;
using namespace std;


template <typename T>
models::ModelState<T,SparseVector> models::SwendsenWang<T>::sample(
	int t,
	models::ModelState<T,SparseVector>& state,
	arithmetic::ThreadOptions& options
) {
	// Compute the temperature and probability of including particular (d-1)-cells.
	double temp = this->parameters.temperatureFunction(t);
	double p = 1-exp(temp);

	// If we're debugging, check that the inclusion probability is actually within
	// acceptable parameters.
	if (this->parameters.DEBUG) {
		cerr << std::format("verifying 0 ≤ {} ≤ 1", p) << endl;
		assert((0 <= p) && (p <= 1));
	}

	// Perform independent Bernoulli trials over each (d-1)-cell to decide whether
	// it gets included in the computation. We do this by performing a (sparse!) matrix
	// multiplication, then checking which of the entries are zero.
	size_t d = this->parameters.dimension;

	SparseVector<T> coefficients = sparse_mat_dot_sparse_vec<T,index_t>(
		this->complex->Coboundary.Matrices[d],
		state.cochain,
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

	// Figure out which indices we're *in*cluding.
	set<int> includeHere;
	
	set_difference(
		this->include.begin(), this->include.end(),
		exclude.begin(), exclude.end(),
		inserter(includeHere, includeHere.begin())
	);

	// Now, sample from the kernel.
	SparseVector<T> sample = arithmetic::submatrixKernelSample<T>(
		this->complex->Coboundary.Matrices[d],	// complete dth coboundary matrix
		this->field,							// field
		exclude,								// rows to exclude
		this->intuniform,						// uniform random over field
		this->RNG,								// RNG
		options,								// compute options
		this->parameters.DEBUG					// debugging
	);

	// If we're debugging, check whether we actually sampled from the kernel.
	if (this->parameters.DEBUG) {
		// Copy and pare down the coboundary matrix.
		SparseMatrix<T> cbd = this->complex->Coboundary.Matrices[d];
		for (auto i : exclude) cbd[i].zero();
		cbd.clear_zero_row();
		cbd.compress();

		// Multiply, and check whether there's anything in the resulting vector;
		// there shouldn't be (i.e. it should have size 0).
		SparseVector<T> outcome = sparse_mat_dot_sparse_vec<T,index_t>(
			cbd,
			sample,
			this->field
		);

		assert(outcome.size() < 1);
	}

	state.cochain = sample;
	state.includes = vector<int>(includeHere.begin(), includeHere.end());
	state.t = t;

	return state;
}

template <typename T>
models::ModelState<T,SparseVector> models::SwendsenWang<T>::initialize(
	models::ModelState<T,SparseVector>& state
) {
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
models::ModelState<T,SparseVector> models::SwendsenWang<T>::initialize(
	SparseVector<T> c,
	models::ModelState<T,SparseVector>& state
) {
	state.cochain = c;
	return state;
}


template <typename T>
models::SwendsenWang<T>::SwendsenWang(
	complexes::Complex<T>* complex,
	ModelParameters parameters
)
	: field(parameters.field > 0 ? Field(SparseRREF::FIELD_Fp, parameters.field) : Field(SparseRREF::FIELD_QQ))
{
	this->parameters = parameters;
	this->complex = complex;

	// Determine the field and build the boundary matrices for the Complex.
	this->complex->constructBoundaryMatrices(this->field);

	// Default set of indices to compare against.
	set<int> _include;
	for (int t=0; t<this->complex->Cells[this->parameters.dimension]; t++) _include.insert(t);

	this->include = _include;

	// Initialize a random number generator.
	std::random_device rd;
	this->RNG = std::mt19937(rd());
	this->unituniform = std::uniform_real_distribution<double>(0,1);
	this->intuniform = std::uniform_int_distribution<int>(0,this->parameters.field > 0 ? this->parameters.field : 1);
}

#endif
