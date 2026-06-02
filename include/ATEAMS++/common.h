
#ifndef ATEAMS_COMMON_H
#define ATEAMS_COMMON_H

#include "libraries/SparseRREF/sparse_mat.h"


#include <set>
#include <vector>
#include <map>

#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <string>


// Utility functions for suppressing output to stderr (specifically because SpaSM
// outputs diagnostic information at every step!).
inline int _suppress() {
	fflush(stderr);
	int fd = dup(STDERR_FILENO);
	freopen("/dev/null", "w", stderr);
	return fd;
}

inline void _resume(int fd) {
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
void printvector(std::vector<t> v) {
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


typedef unsigned long ulong;
typedef uint32_t index_t;
typedef ulong data_t;

typedef SparseRREF::sparse_mat<data_t, index_t> ZpMatrix;
typedef SparseRREF::sparse_vec<data_t, index_t> ZpVector;
typedef SparseRREF::field_t Zp;

typedef std::vector<ZpMatrix> ZpMatrices;
typedef std::vector<ZpVector> ZpVectors;


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
			std::cout << M[i][j] << " ";
		}
		std::cout << std::endl;
	}
}


#endif
