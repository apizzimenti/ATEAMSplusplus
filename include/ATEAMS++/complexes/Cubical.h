
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
		 * 
		 * @tparam RingLike A @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * **Example Usage**
		 * 
		 * @code
		 * using namespace ATEAMS::complexes;
		 * Cubical<Zp> C({3,3,3,3});			// Assumes periodic boundary conditions
		 * Cubical<Zp> D({3,3,3,3}, false); 	// Just a subcomplex of Z^4.
		 * @endcode
		 * 
		 * @var Cubical::name
		 * 	@brief Human-readable short name.
		 * 
		 * @var Cubical::corners
		 * @brief Corners of the lattice. If `corners` is \f$(c_1, c_2, \dots, c_n)\f$,
		 * then the integer lattice created is the product \f$ \prod [0, c_i] \f$; if
		 * periodic boundary conditions are imposed, opposite vertices
		 * (and thus opposite \f$(d-1)\f$-faces) are identified to form
		 * an \f$n\f$-fold torus.
		 */
		template <typename RingLike>
		class Cubical: public Complex<RingLike> {
			public:
				std::string name = "cubical";
				std::vector<int> corners;

				/**
				 * @brief Constructor. See @ref ATEAMS::complexes::Cubical::corners.
				 * 
				 * @param corners Corners of the lattice.
				 * @param periodic Are we imposing periodic boundary conditions?
				 */
				Cubical(std::vector<int> corners, bool periodic);

				/**
				 * @brief Constructor. See @ref ATEAMS::complexes::Cubical::corners.
				 * Imposes periodic boundary conditions by default.
				 * 
				 * @param corners Corners of the lattice.
				 */
				Cubical(std::vector<int> corners);

				void constructBoundaryMatrices(Ring* R) override;
				void constructFlatBoundaryMatrix() override;
				void constructFullBoundaryMatrix(Ring* R) override;
				
				int size() override;
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

#include "ATEAMS++/complexes/Cubical.tpp"

#endif
