
#ifndef ATEAMS_COMPLEXES_CUBICAL_T
#define ATEAMS_COMPLEXES_CUBICAL_T

#ifndef ATEAMS_ARITHMETIC_PERSISTENCE_H
#error __FILE__ should only be included from complexes/Cubical.h.
#endif

#include "ATEAMS++/complexes/Cubical.h"
#include "ATEAMS++/util.h"

#include <cmath>
#include <algorithm>
#include <ranges>
#include <vector>

namespace ATEAMS {
	namespace complexes {
		using namespace std;

		/** @cond */
		template <typename T>
		void _product(vector<vector<T>> sets, vector<vector<T>> &c, vector<T> &active, int depth) {
			for (int i=0; i<sets[depth].size(); i++) {
				active[depth] = sets[depth][i];
				if (depth < sets.size()-1) _product(sets, c, active, depth+1);
				else c.push_back(active);
			}
		}


		template <typename T>
		vector<vector<T>> product(vector<vector<T>> sets) {
			vector<vector<T>> c;
			vector<T> active(sets.size());
			_product<T>(sets, c, active, 0);

			return c;
		}


		template <typename T>
		void _surrounds(vector<T> anchor, vector<vector<T>> &points, vector<T> &active, int depth) {
			for (int f=0; f<2; f++) {
				active[depth] = anchor[depth]+f;
				if (depth < anchor.size()-1) _surrounds(anchor, points, active, depth+1);
				else points.push_back(active);
			}
		}


		template <typename T>
		vector<vector<T>> surrounds(vector<T> anchor) {
			vector<vector<T>> points;
			vector<T> active(anchor.size());
			_surrounds(anchor, points, active, 0);

			return points;
		}


		/**
		 * Creates (the boundary matrices of) a D-dimensional Hamming cube.
		 */
		HammingCube hamming(int D) {
			// All cells in the complex.
			HammingCube complex(D+1);

			// Top-level boundary matrix is just all the vertices.
			vector<bitstring> top(pow(2,D));
			iota(begin(top), end(top), 0);

			vector<vector<bitstring>> cube;
			cube.push_back(top);
			complex[D] = cube;

			// Go from the top dimension down.
			for (int d=D-1; d > -1; d--) {
				// Vector of faces to be filled.
				vector<vector<bitstring>> facesOfDimension;
				vector<vector<bitstring>> cells = complex[d+1];

				/*
				Right now, each "cell" is a list of indices between 0 and (2^D)-1, so
				each d-cell is a list of d indices. Because this is a Hamming cube, we
				can treat these indices as length-D bitstrings; again because this is a
				Hamming cube, the indices of any d-cell are fixed in (D-d) positions.
				For example, in the D=3 Hamming cube, the "front" face has vertices
				[0, 1, 4, 5]:

						(6) 110        111 (7)
							●———————————●
							/|			 /|
					(4)	100	●———————————● |	(5) 101	
							| |	010	(2) | |
							| ●—————————|—● 011 (3)
							|/			|/
							●———————————●
					(0) 000		   001 (1)

				Taking a look at the bitstrings of the vertices composing the front face
				(which is a 2-cell), the second bit is always 0:

							↓
							0 0 0
							0 0 1
							1 0 0
							1 0 1
							↑

				Looking at the "back" face [2, 3, 6, 7], the second bit is always 1.
				Moreover, we can identify the index labels (0, 1, ..., 7) with their 
				_indices_ in the list [0, 1, ..., 7]. Following this pattern, the process
				for identifying the boundary of each d-cell C = [a0, ..., a(d-1)] goes
				like this:

					for each digit p in {0,...,d}:      <--- this is the position of the bit we're flipping
						for each bit b in {0,1}:
							1. create an empty vector F := [] to store the indices of
							the (d-1)-faces of C.
							2. for each index in {0,...,d-1}:
								2.1. set x to be the bitstring i with its pth digit cleared;
									for example, if i = 011 and p = 1, then x := 001.
								2.2. set the pth digit of x to b. if b = 0, this step does
									nothing; if b = 1, then it sets the pth digit to 1.
									for example, if x = 001 and b = 1, then x = 011.
								2.3. add C[x] (i.e. the vertex at index x) to F.
							3. now, the elements of F are the elements of C (i.e. vertex labels)
							at indices whose bitstrings differ in exactly one position
							(i.e. are at Hamming distance 1). because this process flips
							each bit twice, F contains two copies of each label, and
							the labels are unsorted.
							4. remove duplicates and sort the elements of F. these are
							the (d-1)-dimensional faces of C.


				For example, let C = [0, 1, 4, 5]. The vertices are in lexicographic order,
				so we know that labels at _indices_ differing in (d-1) bits belong to the
				same face. (If d=3 --- i.e. we're looking for the 2-faces of the cube --- then
				the vertex labels whose indices share (3-2)=1 bit and differ in d=2 bits
				belong to the same face.) Applying our indexing logic, we get the following
				index pairs indicating edges:

								00 01 ———> [ C[0], C[1] ] ———> [ 0, 1 ]
								00 10 ———> [ C[0], C[2] ] ———> [ 0, 4 ]
								01 11 ———> [ C[1], C[3] ] ———> [ 1, 4 ]
								10 11 ———> [ C[2], C[3] ] ———> [ 4, 5 ]

				which are exactly the edges bounding the front face of C. Thus, the boundary
				(in terms of vertex labels) of C is

								[ [ 0, 1 ], [ 0, 4 ], [ 1, 4 ], [ 4, 5 ] ].

				(In practice, each of these edges shows up twice in F, so we just sort F
				and then delete every other entry.)
				*/
				for (auto cell : cells) {
					for (int pos=0; pos < d+1; pos++) {
						for (char flip=0; flip<2; flip++) {
							vector<bitstring> faces;
							for (bitstring i=0; i < cell.size(); i++) {
								bitstring x = i & ~(1 << pos);
								x = x | (flip << pos);
								faces.push_back(cell[x]);
							}

							// Sort and take every other entry.
							sort(faces.begin(), faces.end());
							vector<bitstring> face;
							for (int i=0; i<faces.size(); i+=2) face.push_back(faces[i]);

							facesOfDimension.push_back(face);
						}
					}
				}

				// Remove duplicates; lexicographic sort.
				sort(facesOfDimension.begin(), facesOfDimension.end());
				auto it = unique(facesOfDimension.begin(), facesOfDimension.end());
				facesOfDimension.erase(it, facesOfDimension.end());

				complex[d] = facesOfDimension;
			}

			return complex;
		}


		/**
		 * Creates flat boundary matrices for a HammingCube.
		 */
		HammingCubeBoundary hammingBoundary(HammingCube cube) {
			int D = cube.size();
			HammingCubeBoundary Boundary(D);

			for (int d=D-1; d > 0; d--) {
				vector<vector<bitstring>> cells = cube[d];
				vector<vector<bitstring>> faces = cube[d-1];
				vector<vector<int>> boundaries;

				for (auto cell : cells) {
					vector<int> boundary;
					for (int i=0; i<faces.size(); i++) {
						bool hasFace = includes(cell.begin(), cell.end(), faces[i].begin(), faces[i].end());
						if (hasFace) boundary.push_back(i);
					}
					boundaries.push_back(boundary);
				}
				Boundary[d] = boundaries;
			}

			vector<int> _zero(1);
			vector<vector<int>> zero(cube[0].size(),_zero);
			Boundary[0] = zero;

			return Boundary;
		}

		/**
		 * Creates a vector of indexers --- std::maps that take vectors to indices --- to
		 * construct the full boundary matrix.
		 */
		vector<indexer> protoCubicalLattice(vector<int> corners, HammingCube Cube, HammingCubeBoundary Boundary, bool periodic) {
			vector<indexer> proto(corners.size()+1);

			// Get the "axes" of the space.
			vector<vector<int>> axes;
			for (auto l : corners) {
				vector<int> axis(l);
				iota(begin(axis), end(axis), 0);
				axes.push_back(axis);
			}

			// Get all the corner points.
			vector<vector<int>> coordinates = ATEAMS::complexes::product<int>(axes);

			// For each corner point, figure out its surrounding points, then determine
			// faces/edges/etc. If we're using periodic boundary conditions, we need to
			// re-set any coordinate that contains a boundary point to its corresponding
			// opposite point.
			for (int i=0; i < coordinates.size(); i++) {
				HammingCubeBoundary localBoundary(proto.size());

				// Get the anchor (bottom-left) point of the cube, and find the indices
				// of the points surrounding it in lexicographic order.
				vector<int> anchor = coordinates[i];
				vector<vector<int>> localVertices = surrounds<int>(anchor);
				vector<int> localIndices(localVertices.size());

				for (int j=0; j<localVertices.size(); j++) {
					// With periodic boundary conditions, scan the vertices to check
					// whether it's an extremal point (i.e. one of the coordinates is a
					// corner coordinate).
					if (periodic) {
						for (int k=0; k<localVertices[j].size(); k++) {
							if (localVertices[j][k] == corners[k]) localVertices[j][k] = 0;
						}
					}

					// Insert the points if they haven't been inserted yet.
					proto[0].try_emplace(localVertices[j], proto[0].size());
					localIndices[j] = proto[0][localVertices[j]];
				}

				localBoundary[0] = localVertices;
				
				// From here, repeat a similar process, working up in dimension.
				for (int d=1; d<proto.size(); d++) {
					vector<vector<int>> canonicalBoundary = Boundary[d];
					vector<vector<int>> localCellBoundary;

					for (auto b : canonicalBoundary) {
						vector<int> localCell(b.size());
						for (int j=0; j<b.size(); j++) localCell[j] = proto[d-1][localBoundary[d-1][b[j]]];

						proto[d].try_emplace(localCell, proto[d].size());
						localCellBoundary.push_back(localCell);
					}

					localBoundary[d] = localCellBoundary;
				}
			}

			return proto;
		}

		/**
		 * From a proto-boundary (a vector of indexers), construct the full boundary.
		 */
		Lattice cubicalLattice(vector<indexer> proto) {
			Lattice cubical(proto.size());

			// Each of the boundaries should be in the proto-boundary, with an
			// assigned index.
			for (int d=1; d<proto.size(); d++) {
				int count = proto[d].size();
				cubical[d] = vector<vector<int>>(count);

				for (auto &[b, i] : proto[d]) {
					cubical[d][i] = b;
				}
			}

			// Fill in the vertices, since they're otherwise indexed by coordinates.
			vector<vector<int>> vertices(proto[0].size());
			for (int i=0; i<vertices.size(); i++) vertices[i] = {i};
			cubical[0] = vertices;

			return cubical;
		}


		/**
		 * Using the full (flat) boundary matrices, make SparseRREF matrices.
		 */
		template <typename T>
		SparseMatrices<T> sparseBoundaryMatrices(Lattice L, Field F) {
			SparseMatrices<T> Boundary(L.size());
			int mod = (int)F.mod.n;

			for (int d=L.size()-1; d > 0; d--) {
				// Get the dimensions of the matrix; instantiate.
				int M = L[d-1].size(), N = L[d].size();
				SparseMatrix<T> B(M,N);

				// Fill in the values.
				for (int col=0; col < L[d].size(); col++) {
					for (int i=0; i < L[d][col].size(); i++) {
						int row = L[d][col][i];
						SparseVector<T>& matrow = B.rows[row];
						matrow.push_back(
							(index_t)col,
							(T)((bool)(i%2) ? mod-1 : 1)
						);
					}
				}

				// Push the matrix to the collection.
				B.compress();
				Boundary[d] = B;
			}

			return Boundary;
		}


		FlatBoundaryMatrix flatBoundaryMatrix(Lattice L, vector<int> offsets) {
			// Construct the full flat boundary matrix.
			FlatBoundaryMatrix flat;
			int offset;

			for (int d=0; d < L.size(); d++) {
				offset = (d > 1) ? offsets[d-2] : 0;

				for (int c=0; c < L[d].size(); c++) { // haha c++
					vector<int> cell;

					if (d > 0) {
						cell = L[d][c];
						for (int j=0; j < cell.size(); j++) cell[j] = cell[j] + offset;
					}
					
					flat.push_back(cell);
				}
			}

			return flat;
		}


		template <typename T>
		SparseMatrix<T> sparseFullBoundaryMatrix(FlatBoundaryMatrix flat, Field F) {
			SparseMatrix<T> Full(flat.size(), flat.size());
			int mod = (int)F.mod.n;

			vector<int> cell;

			for (int col=0; col < flat.size(); col++) {
				cell = flat[col];
				for (int i=0; i < cell.size(); i++) {
					int row = flat[col][i];
					SparseVector<T>& matrow = Full.rows[row];
					matrow.push_back(
						(index_t)col,
						(T)((bool)(i%2) ? mod-1 : 1)
					);
				}
			}

			Full.compress();
			return Full;
		}


		vector<int> getCellCounts(Lattice L) {
			vector<int> counts;
			for (int d=0; d < L.size(); d++) counts.push_back(L[d].size());

			return counts;
		}

		vector<int> makeOffsets(vector<int> Cells) {
			vector<int> offsets(Cells.size());
			std::partial_sum(Cells.begin(), Cells.end(), offsets.begin());

			return offsets;
		}

		vector<vector<int>> makeBreaks(vector<int> counts) {
			// Get cell counts and breaks.
			vector<vector<int>> breaks(counts.size());
			int last = 0;

			for (int d=0; d < counts.size(); d++) {
				vector<int> broken(2);
				broken[0] = last;
				broken[1] = broken[0] + counts[d];
				breaks[d] = broken;

				last = broken[1];
			}

			return breaks;
		}


		/** @endcond */


		template <typename T>
		int Cubical<T>::size() {
			return this->_size;
		}


		template <typename T>
		void Cubical<T>::constructFlatBoundaryMatrix() {
			// Create a Hamming cube of the appropriate dimension.
			HammingCube cube = hamming(this->corners.size());
			HammingCubeBoundary cubeBoundary = hammingBoundary(cube);

			// Compute the proto-boundary (which includes a vertex map), and get the
			// flat boundary matrices for each dimension.
			vector<indexer> protoBoundary = protoCubicalLattice(corners, cube, cubeBoundary, this->periodic);
			Lattice L = cubicalLattice(protoBoundary);

			// Get cell counts and breaks.
			this->Cells = getCellCounts(L);
			this->_size = std::accumulate(this->Cells.begin(), this->Cells.end(), 0);
			this->breaks = makeBreaks(this->Cells);
			this->offsets = makeOffsets(this->Cells);

			this->Boundary.Flat = flatBoundaryMatrix(L, this->offsets);
		}


		template <typename T>
		void Cubical<T>::constructBoundaryMatrices(Field F) {
			// Create a Hamming cube of the appropriate dimension.
			HammingCube cube = hamming(this->corners.size());
			HammingCubeBoundary cubeBoundary = hammingBoundary(cube);

			// Compute the proto-boundary (which includes a vertex map), and get the
			// flat boundary matrices for each dimension.
			vector<indexer> protoBoundary = protoCubicalLattice(corners, cube, cubeBoundary, this->periodic);
			Lattice L = cubicalLattice(protoBoundary);

			// Get cell counts and breaks.
			this->Cells = getCellCounts(L);
			this->_size = std::accumulate(this->Cells.begin(), this->Cells.end(), 0);
			this->breaks = makeBreaks(this->Cells);
			this->offsets = makeOffsets(this->Cells);

			// Get the sparse (co)boundary matrices.
			this->Boundary.Matrices = sparseBoundaryMatrices<T>(L, F);
			this->Coboundary.Matrices = SparseMatrices<T>(this->Boundary.Matrices.size());

			for (int d=0; d < this->Boundary.Matrices.size(); d++) {
				this->Coboundary.Matrices[d] = this->Boundary.Matrices[d].transpose();
				this->Coboundary.Matrices[d].compress();
			}
		}


		template <typename T>
		void Cubical<T>::constructFullBoundaryMatrix(Field F) {
			// Create a Hamming cube of the appropriate dimension.
			HammingCube cube = hamming(this->corners.size());
			HammingCubeBoundary cubeBoundary = hammingBoundary(cube);

			// Compute the proto-boundary (which includes a vertex map), and get the
			// flat boundary matrices for each dimension.
			vector<indexer> protoBoundary = protoCubicalLattice(corners, cube, cubeBoundary, this->periodic);
			Lattice L = cubicalLattice(protoBoundary);

			// Get cell counts and breaks.
			this->Cells = getCellCounts(L);
			this->_size = std::accumulate(this->Cells.begin(), this->Cells.end(), 0);
			this->breaks = makeBreaks(this->Cells);
			this->offsets = makeOffsets(this->Cells);

			FlatBoundaryMatrix flat = flatBoundaryMatrix(L, this->offsets);
			this->Boundary.Full = sparseFullBoundaryMatrix<T>(flat, F);
			this->Coboundary.Full = this->Boundary.Full.transpose();
		}


		template <typename T>
		Cubical<T>::Cubical(vector<int> corners, bool periodic) {
			this->corners = corners;
			this->periodic = periodic;
		}

		template <typename T>
		Cubical<T>::Cubical(vector<int> corners) {
			this->corners = corners;
			this->periodic = true;
		}
	}
}

#endif

