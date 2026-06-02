
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/SwendsenWang.h"
#include "libraries/SparseRREF/sparse_mat.h"

#include <random>
#include <algorithm>

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

	// The zero entries are those that do *not* show up in the array of .indices.
	// TODO it's unclear whether this is efficient at all!! Probably isn't, so
	// should change.
	int N = this->complex->Cells[d];
	vector<int> include(N);
	iota(begin(include), end(include), 0);

	if (coefficients.size() > 0) {
		for (int i=0; i < coefficients.size(); i++) {
			include.erase(remove(include.begin(), include.end(), coefficients.indices[i]), include.end());
		}
	}

	// Now, get a submatrix view

	ZpVector r;

	return r;
}

/**
 * Creates an initial cochain of all zeros.
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
	this->RNG = std::mt19937_64(rd());
	this->unituniform = std::uniform_real_distribution<double>(0,1);
	this->intuniform = std::uniform_int_distribution<int>(0,this->parameters.field);
}
