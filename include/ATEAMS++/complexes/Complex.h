
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
		 * @tparam RingLike A @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @var BoundaryData::Matrices
		 * 	@ref SparseMatrices where `Matrices[d]` is the \f$d\f$th boundary
		 * 	matrix \f$\partial_d\f$ of the given complex.
		 * 
		 * @var BoundaryData::Bases
		 * 	@ref SparseBases where `Bases[d]` is a @ref SparseBasis representing
		 * 	a basis for the \f$d\f$th homology group \f$H_d(X)\f$.
		 * 
		 * @var BoundaryData::Full
		 * 	A @ref SparseMatrix representing the full boundary matrix \f$\partial\f$
		 * 	for the complex \f$X\f$ --- that is, the square matrix with rows and
		 * 	columns indexed by all cells of \f$X\f$ where \f$\partial_{ij} = \pm 1\f$
		 * 	if cell \f$i\f$ is a face of cell \f$j\f$, and \f$0\f$ otherwise.
		 * 
		 * @var BoundaryData::Flat
		 * 	A vector of vectors representing the fully indexed flat boundary
		 * 	matrix. For example, the \f$d\f$th entry contains the indices of its
		 * 	\f$(d-1)\f$-dimensional faces.
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
		 * @brief Template class for various complexes.
		 * @tparam RingLike A @ref Ring, like @ref Zp or @ref Q.
		 * 
		 * @var Complex::Boundary
		 * 	@brief Boundary data.
		 * 
		 * @var Complex::Coboundary
		 * 	@brief Coboundary data.
		 * 
		 * @var Complex::Cells
		 * 	@brief Cell counts; `Cells[d]` is the number of \f$d\f$-dimensional cells
		 * 	in the complex.
		 * 
		 * @var Complex::Breaks
		 * 	@brief The \f$d\f$th entry indicates the starting and ending
		 *  indices of \f$d\f$-cells in the full/flat boundary matrices.
		 * 
		 * @var Complex::Offsets
		 * 	@brief The \f$d\f$th entry indicates the cumulative number of
		 *  cells over all dimensions lower than d.
		 * 
		 * @var Complex::_size
		 * 	@brief Number of cells in the complex.
		 * 
		 * @var Complex::periodic
		 * 	@brief Does this complex have periodic boundary conditions?
		 */
		template <typename RingLike>
		class Complex {
			public:
				BoundaryData<RingLike> Boundary;
				BoundaryData<RingLike> Coboundary;
				
				std::vector<int> Cells;
				std::vector<std::vector<int>> Breaks;
				std::vector<int> Offsets;

				int _size;
				bool periodic;

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
				 * @param R (Pointer to a) coefficient ring, like @ref Zp or @ref Q.
				 */
				virtual void constructBoundaryMatrices(Ring* R) = 0;

				/**
				 * @brief Constructs flat boundary matrix.
				 */
				virtual void constructFlatBoundaryMatrix() = 0;

				/**
				 * @brief Constructs a full boundary matrix (in SparseRREF format).
				 * @param R (Pointer to a) coefficient ring, like @ref Zp or @ref Q.
				 */
				virtual void constructFullBoundaryMatrix(Ring* R) = 0;

				/** @brief Total number of cells. */
				virtual int size() = 0;
		};
	}
}

#endif
