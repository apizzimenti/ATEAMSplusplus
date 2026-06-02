
#ifndef ATEAMS_COMPLEXES_CUBICAL_H
#define ATEAMS_COMPLEXES_CUBICAL_H

#include "ATEAMS++/complexes/Complex.h"

#include <vector>
#include <string>
#include <map>

namespace ATEAMS {
	class Cubical: public Complex {
		public:
			std::string kind = "cubical";							// what kind of complex is this?
			std::vector<int> corners;								// "corners" of the lattice.

			Cubical(std::vector<int> corners, bool periodic);		// constructor

			void constructBoundaryMatrices(Zp F) override;		// construct boundary matrices.
	};


	// Type declarations used only for constructing cubical complexes.
	typedef uint16_t bitstring;
	typedef std::vector<std::vector<std::vector<bitstring>>> HammingCube;

	typedef std::vector<std::vector<std::vector<int>>> HammingCubeBoundary;
	typedef HammingCubeBoundary Lattice;
	typedef std::map<std::vector<int>,int> indexer;
}

#endif
