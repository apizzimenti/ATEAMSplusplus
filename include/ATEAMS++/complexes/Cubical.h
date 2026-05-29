
#ifndef ATEAMS_COMPLEX_CUBICAL_H
#define ATEAMS_COMPLEX_CUBICAL_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/complexes/Complex.h"

class Cubical: public Complex {
	public:
		string kind = "cubical";				// what kind of complex is this?

		vector<int> corners;					// "corners" of the lattice.

		Cubical(vector<int> corners, bool periodic);
};

#endif
