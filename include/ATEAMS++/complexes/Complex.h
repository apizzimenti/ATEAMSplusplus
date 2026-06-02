
#ifndef ATEAMS_COMPLEXES_COMPLEX_H
#define ATEAMS_COMPLEXES_COMPLEX_H

#include "ATEAMS++/common.h"

#include <vector>
#include <string>


namespace ATEAMS {
	// Boundary data.
	struct BoundaryData {
		ZpMatrices Matrices;		// vector of sparse boundary matrices;
		ZpVectors Basis;			// vector of sparse vectors (cycle basis);
		ZpMatrix Full;				// flattened full boundary matrix.
		std::vector<int> breaks;	// indices of .Full where dimensions change.
	};

	class Complex {
		public:
			// IO
			static void toFile(std::string filename);
			static void fromFile(std::string filename);

			// Build boundary matrices.
			virtual void constructBoundaryMatrices(Zp F) = 0;

			BoundaryData Boundary;			// (Co)boundary information.
			BoundaryData Coboundary;

			// Diagnostic info?
			std::vector<int> Cells;				// counts cells per dimension
			bool periodic;					// is this complex periodic?

			virtual ~Complex() = default;
	};
}

#endif
