
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
		 * 	@ref SparseMatrices where `Matrices[d]` is the \f$d\f$th boundary
		 * 	matrix \f$\partial_d\f$ of the given complex.
		 * 
		 * @var BoundaryData::Bases
		 * @ref SparseBases where `Bases[d]` is a @ref SparseBasis representing
		 * 	a basis for the \f$d\f$th homology group \f$H_d(X)\f$.
		 * 
		 * @var BoundaryData::Full
		 * A @ref SparseMatrix representing the full boundary matrix \f$\partial\f$
		 * 	for the complex \f$X\f$ --- that is, the square matrix with rows and
		 * 	columns indexed by all cells of \f$X\f$ where \f$\partial_{ij} = \pm 1\f$
		 * 	if cell \f$i\f$ is a face of cell \f$j\f$, and \f$0\f$ otherwise.
		 * 
		 * @var BoundaryData::Flat
		 * A vector of vectors representing the fully indexed flat boundary
		 * matrix. For example, the \f$d\f$th entry contains the indices of its
		 * \f$(d-1)\f$-dimensional faces.
		 */
		template <typename RingLike>
		struct BoundaryData {
			SparseMatrices<RingLike> Matrices;
			SparseBases<RingLike> Bases;
			SparseMatrix<RingLike> Full;
			FlatBoundaryMatrix Flat;
		};

		/**
		 * @class Complex
		 * @brief Abstract (template) class for various complexes.
		 */
		template <typename RingLike>
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
				 * @param F Finite field \f$ \Z/p\Z \f$ or \f$ \Q \f$..
				 */
				virtual void constructBoundaryMatrices(Ring* R) = 0;

				/**
				 * @brief Constructs flat boundary matrix.
				 */
				virtual void constructFlatBoundaryMatrix() = 0;

				/**
				 * @brief Constructs a full boundary matrix (in SparseRREF format).
				 */
				virtual void constructFullBoundaryMatrix(Ring* R) = 0;

				/** @brief Total number of cells. */
				virtual int size() = 0;

				/** @brief Boundary data. */
				BoundaryData<RingLike> Boundary;

				/** @brief Coboundary data. */
				BoundaryData<RingLike> Coboundary;

				/** @brief Contains cell counts. */
				std::vector<int> Cells;

				/**
				 * @brief The \f$d\f$th entry indicates the starting and ending
				 * indices of \f$d\f$-cells in the full/flat boundary matrices.
				 */
				std::vector<std::vector<int>> breaks;

				/**
				 * @brief The \f$d\f$th entry indicates the cumulative number of
				 * cells over all dimensions lower than d.
				 */
				std::vector<int> offsets;

				/** @brief Total cell count. */
				int _size;

				/** @brief Periodic boundary conditions? */
				bool periodic;
		};
	}
}

#endif
