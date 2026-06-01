
#ifndef ATEAMS_COMMON_H
#define ATEAMS_COMMON_H

class ATEAMS {
	public:
		void print();
};


#include <set>
#include <vector>
#include <map>

#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <string>

using namespace std;


// Typedefs for SparseRREF (used in Sampling and Persistence).
typedef unsigned long ulong;
typedef uint32_t index_t;
typedef ulong data_t;

// Type for indices and type for data; for now, we use unsigned 32-bit ints for
// indices (since we're never negative) and signed chars for data (since our
// fields are never that big).
typedef int32_t INDEXTYPE;
typedef char DATATYPE;

// Standard sets, indices
typedef set<INDEXTYPE> Set;
typedef map<INDEXTYPE,INDEXTYPE> Map;
typedef vector<INDEXTYPE> Index;

// Arithmetic.
typedef vector<DATATYPE> Lookup;
typedef vector<Lookup> Table;

// Persistence pairs; these are literally just pairs of ints.
typedef vector<Index> PersistencePairs;

// Matrices.
typedef map<INDEXTYPE,DATATYPE> Column;
typedef vector<Column> BoundaryMatrix;
typedef vector<Index> FlatBoundaryMatrix;

typedef vector<Column> Basis;
typedef vector<Basis> Bases;

// Utility functions for suppressing output to stderr (specifically because SpaSM
// outputs diagnostic information at every step!).
int _suppress() {
	fflush(stderr);
	int fd = dup(STDERR_FILENO);
	freopen("/dev/null", "w", stderr);
	return fd;
}

void _resume(int fd) {
	fflush(stderr);
	dup2(fd, fileno(stderr));
	close(fd);
}


// Helper method for printing `map`s.'
template <typename MapStorage>
void printmap(MapStorage m) {
	std::cout << "{ ";
	for (const auto& [k, v] : m) {
		std::cout << (int)k << ": " << (int)v << ", ";
	}
	std::cout << " }" << std::endl;
}

// Helper method for printing `vector`s.
template <typename t>
void printvector(vector<t> v) {
	std::cout << "[ ";
	for (auto &k : v) {
		std::cout << k << " ";
	}
	std::cout << "]" << std::endl;
}


template <typename MapStorage>
void printvectormap(MapStorage m) {
	std::cout << "{ ";
	for (const auto& [k, v] : m) {
		std::cout << "[ ";
		for (auto &p : k) {
			std::cout << p << " ";
		}
		std::cout << "]";
		
		std::cout << ": " << (int)v << ", ";
	}
	std::cout << " }" << std::endl;
}


#include "libraries/SparseRREF/sparse_mat.h"

typedef SparseRREF::sparse_mat<data_t, index_t> ZpMatrix;
typedef SparseRREF::sparse_vec<data_t, index_t> ZpVector;
typedef SparseRREF::field_t Zp;

typedef vector<ZpMatrix> ZpMatrices;
typedef vector<ZpVector> ZpVectors;


// Helper method for printing SparseRREF matrices.
template <typename Matrix>
void printSparseRREFmat(Matrix A) {
	int rows, columns;
	rows = A.nrow;
	columns = A.ncol;

	std::vector<std::vector<int>> M(rows, std::vector<int>(columns, 0));

	for (int i=0; i<rows; i++) {
		for (auto [j,v] : A.rows[i]) {
			M[i][j] = (int)v;
		}
	}

	for (int i=0; i<rows; i++) {
		for (int j=0; j<columns; j++) {
			cout << M[i][j] << " ";
		}
		cout << endl;
	}
}


#endif
