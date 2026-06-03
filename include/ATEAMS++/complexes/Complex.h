
#ifndef ATEAMS_COMPLEXES_COMPLEX_H
#define ATEAMS_COMPLEXES_COMPLEX_H

#include "ATEAMS++/common.h"

#include <vector>
#include <string>


namespace ATEAMS {
	namespace complexes {
		/**
		 * @struct BoundaryData
		 * @brief Encodes a complex by its boundary matrices.
		 * 
		 * @var BoundaryData::Matrices
		 * @ref ATEAMS::ZpMatrices where `Matrices[d]` is the \f$d\f$th boundary
		 * 	matrix \f$\partial_d\f$ of the given complex.
		 * 
		 * @var BoundaryData::Bases
		 * @ref ATEAMS::Bases where `Bases[d]` is a @ref ATEAMS::ZpBasis representing
		 * 	a basis for the \f$d\f$th homology group \f$H_d(X)\f$.
		 * 
		 * @var BoundaryData::Full
		 * A @ref ATEAMS::ZpMatrix representing the full boundary matrix \f$\partial\f$
		 * 	for the complex \f$X\f$ --- that is, the square matrix with rows and
		 * 	columns indexed by all cells of \f$X\f$ where \f$\partial_{ij} = \pm 1\f$
		 * 	if cell \f$i\f$ is a face of cell \f$j\f$, and \f$0\f$ otherwise.
		 * 
		 * @var BoundaryData::breaks
		 * The \f$d\f$th entry indicates the starting index of \f$d\f$-cells in
		 * 	the full boundary matrix.
		 */
		struct BoundaryData {
			ZpMatrices Matrices;
			ZpBases Bases;
			ZpMatrix Full;
			std::vector<int> breaks;
		};

		/**
		 * @class Complex
		 * @brief Abstract (template) class for various complexes.
		 */
		class Complex {
			public:
				/**
				 * @brief Writes boundary matrices to file.
				 * @param filename Destination.
				 */
				static void toFile(std::string filename);

				/**
				 * @brief Reads boundary matrices from file.
				 * @param filename Source.
				 */
				static void fromFile(std::string filename);

				/**
				 * @brief Constructs boundary matrices.
				 * @param F Finite field \f$ \Z/p\Z \f$.
				 */
				virtual void constructBoundaryMatrices(Zp F) = 0;

				/** @brief Boundary data. */
				BoundaryData Boundary;

				/** @brief Coboundary data. */
				BoundaryData Coboundary;

				/** @brief Contains cell counts. */
				std::vector<int> Cells;

				/** @brief Periodic boundary conditions? */
				bool periodic;
		};
	}
}

#endif
