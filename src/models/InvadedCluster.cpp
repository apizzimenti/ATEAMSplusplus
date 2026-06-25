
#include "ATEAMS++/util.h"
#include "ATEAMS++/common.h"
#include "ATEAMS++/models/InvadedCluster.h"
#include "ATEAMS++/arithmetic/kernel.h"
#include "ATEAMS++/arithmetic/persistence.h"

#include <SparseRREF/sparse_vec.h>

#include <random>
#include <algorithm>
#include <cassert>

using namespace ATEAMS;
using namespace std;


ZpVector models::InvadedCluster::sample(int t, arithmetic::ThreadOptions& options) {
	int d = this->parameters.dimension;

	// Determine which (d-1)-cells are "satisfiable" (i.e. those on which the
	// current (d-1)-cochain vanishes).
	ZpVector coefficients = sparse_mat_dot_sparse_vec<data_t, index_t>(
		this->complex->Coboundary.Matrices[d],
		this->state.cochain,
		this->field
	);

	// The zero entries are those that don't show up in the indices, so we
	// automatically exclude these.
	int N = this->complex->Cells[d];
	set<size_t> exclude;
	vector<size_t> include;
	
	if (coefficients.size() > 0) {
		for (int e=0; e < coefficients.size(); e++) exclude.insert((size_t)e);
	}

	for (int i=0; i < N; i++) {
		if (!exclude.contains(i)) include.push_back(i);
	}

	std::shuffle(include.begin(), include.end(), this->RNG);

	// Construct the filtration.
	vector<size_t> inside(include.begin(), include.end());
	inside.insert(inside.end(), exclude.begin(), exclude.end());

	int stop = this->complex->breaks[d][1];
	int start = this->complex->breaks[d][0];
	
	// Fill in all the indices for cells of dimension not equal to d.
	// TODO a small optimization: we only have to do the persistence algorithm
	// through the (d+1)-dimensional cells to capture the d-dimensional events.
	// This will probably save a boatload of time!
	for (int j=0; j < inside.size(); j++) this->filtration[j+start] = inside[j]+start;

	// Persist.
	vector<int> essential;

	if (this->parameters.field < 3) {
		essential = arithmetic::PHATPersistence(this->complex, this->filtration, this->parameters.dimension);
	} else {
		essential = arithmetic::TwistPersistence(this->complex, this->filtration, this->field, this->parameters.dimension);
	}

	std::erase_if(essential, [stop, start](int x) { return !((start <= x) && (x < stop)); });
	std::sort(essential.begin(), essential.end());

	if (this->parameters.DEBUG) {
		cerr << "filtration:" << endl;
		printvector<int>(this->filtration);
		cerr << endl;

		cerr << std::format("essential cycles (of dimension) {}", d) << endl;
		printvector<int>(essential);
		cerr << endl;
	}

	// Now, sample from the kernel, excluding unoccupied columns *or* columns
	// corresponding to cells included *after* percolating.
	int stopat = essential[this->parameters.stoppingFunction(t)];
	
	set<size_t> leaveout;
	vector<int> included(stopat-start, 0);

	for (int i=stopat; i < stop; i++) leaveout.insert(this->filtration[i]-start);
	for (int i=start; i < stopat; i++) included[i-start] = this->filtration[i]-start;

	ZpVector sample = arithmetic::submatrixKernelSample(
		this->complex->Coboundary.Matrices[d],	// complete dth coboundary matrix
		this->field,							// field
		leaveout,								// rows to exclude
		this->intuniform,						// uniform random over field
		this->RNG,								// RNG
		options,								// compute options
		this->parameters.DEBUG					// debugging
	);

	// If we're debugging, check whether we actually sampled from the kernel.
	if (this->parameters.DEBUG) {
		// Copy and pare down the coboundary matrix.
		ZpMatrix cbd = this->complex->Coboundary.Matrices[d];
		for (auto i : leaveout) cbd[i].zero();
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
	this->state.includes = included;

	return sample;
}


void models::InvadedCluster::initialize() {
	size_t dimension = this->parameters.dimension-1;
	int N = this->complex->Cells[dimension];

	ZpVector cochain;
	for (int i=0; i < N; i++) cochain.push_back(
		(index_t)i, (data_t)this->intuniform(this->RNG)
	);

	cochain.compress();
	this->state.cochain = cochain;
}


void models::InvadedCluster::initialize(ZpVector c) {
	this->state.cochain = c;
}

models::InvadedCluster::InvadedCluster(complexes::Complex* complex, InvadedClusterParameters parameters)
	: field(Zp(SparseRREF::FIELD_Fp, parameters.field))
{
	this->parameters = parameters;
	this->complex = complex;

	// Determine the field and build the boundary matrices for the Complex.
	this->complex->constructBoundaryMatrices(this->field);

	if (this->parameters.field < 3) this->complex->constructFlatBoundaryMatrix();
	else this->complex->constructFullBoundaryMatrix(this->field);

	vector<int> filtration(this->complex->size());
	std::iota(filtration.begin(), filtration.end(), 0);
	this->filtration = filtration;

	// Initialize a random number generator.
	std::random_device rd;
	this->RNG = std::mt19937(rd());
	this->intuniform = std::uniform_int_distribution<int>(0,this->parameters.field);
}
