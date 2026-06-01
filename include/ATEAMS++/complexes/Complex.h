
#ifndef ATEAMS_COMPLEX_H
#define ATEAMS_COMPLEX_H

#include "ATEAMS++/common.h"

using namespace std;


class Complex {
	public:
		// IO
		static void toFile(string filename);
		static void fromFile(string filename);

		// Boundary data.
		struct BoundaryData {
			ZpMatrices Matrices;		// vector of sparse boundary matrices;
			ZpVectors Basis;			// vector of sparse vectors (cycle basis);
			ZpMatrix Full;				// flattened full boundary matrix.
			vector<int> breaks;			// indices of .Full where dimensions change.
		};

		BoundaryData Boundary;

		// Coboundary data.
		struct CoboundaryData {
			ZpMatrices Matrices;		// vector of sparse coboundary matrices;
			ZpVectors Basis;			// vector of sparse covectors (cocycle basis);
			ZpMatrix Full;				// flattened full coboundary matrix.
			vector<int> breaks;			// indices of .Full where dimensions change.
		};

		CoboundaryData Coboundary;

		// Diagnostic info?
		vector<int> Cells;				// counts cells per dimension
		bool periodic;					// is this complex periodic?
		
		// Field, since we need this to build the matrices.
		Zp F;
};

#endif
