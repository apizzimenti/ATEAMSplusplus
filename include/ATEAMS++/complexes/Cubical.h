
#ifndef ATEAMS_COMPLEXES_CUBICAL_H
#define ATEAMS_COMPLEXES_CUBICAL_H

#include "ATEAMS++/common.h"
#include "ATEAMS++/complexes/Complex.h"

#include <vector>
#include <string>
#include <map>

namespace ATEAMS {
	namespace complexes {
		/**
		 * @class Cubical
		 * @brief Cubical complex.
		 */
		class Cubical: public Complex {
			public:
				/** @brief Short name, mostly for metadata collection. */
				std::string kind = "cubical";

				/**
				 * @brief Corners of the lattice.
				 * 
				 * If `corners` is \f$(c_1, c_2, \dots, c_n)\f$, then the integer
				 * lattice created is the product \f$ \prod [0, c_i] \f$; if
				 * periodic boundary conditions are imposed, opposite vertices
				 * (and thus opposite \f$(d-1)\f$-faces) are identified to form
				 * an \f$n\f$-fold torus.
				 */
				std::vector<int> corners;

				/**
				 * @brief Constructor. See @ref ATEAMS::complexes::Cubical::corners.
				 * 
				 * @param corners Corners of the lattice.
				 * @param periodic Are we imposing periodic boundary conditions?
				 */
				Cubical(std::vector<int> corners, bool periodic);

				/**
				 * @brief Constructs boundary matrices for the cubical complex
				 * with specified corners and boundary conditions.
				 */
				void constructBoundaryMatrices(Zp F) override;
		};

		/** @cond */
		typedef uint16_t bitstring;
		typedef std::vector<std::vector<std::vector<bitstring>>> HammingCube;

		typedef std::vector<std::vector<std::vector<int>>> HammingCubeBoundary;
		typedef HammingCubeBoundary Lattice;
		typedef std::map<std::vector<int>,int> indexer;
		/** @endcond */
	}
}

#endif
