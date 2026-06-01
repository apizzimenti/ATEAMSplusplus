
#include "ATEAMS++/complexes/Cubical.h"
#include "ATEAMS++/util.h"

#include <cmath>
#include <algorithm>
#include <ranges>

using namespace std;


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


typedef map<vector<int>,int> indexer;

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
	vector<vector<int>> coordinates = product<int>(axes);

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
			// whether it's a "boundary" point.
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

		vector<vector<int>> localEdges;

		// Now, given local indices, we can create edges and insert into the
		// prototype boundary.
		for (auto edge : Boundary[1]) {
			vector<int> localEdge(2);
			for (int j=0; j<edge.size(); j++) localEdge[j] = proto[0][localVertices[edge[j]]];

			proto[1].try_emplace(localEdge, proto[1].size());
			localEdges.push_back(localEdge);
		}

		localBoundary[1] = localEdges;

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


Cubical::Cubical(vector<int> corners, bool periodic=true) {
	this->corners = corners;
	this->periodic = periodic;

	// Create a Hamming cube of the appropriate dimension.
	HammingCube cube = hamming(corners.size());
	HammingCubeBoundary cubeBoundary = hammingBoundary(cube);

	// Compute the proto-boundary (which includes a vertex map), and get the
	// flat boundary matrices for each dimension.
	vector<indexer> protoBoundary = protoCubicalLattice(corners, cube, cubeBoundary, periodic);
	Lattice Boundary = cubicalLattice(protoBoundary);

	// Get cell counts.
	for (int d=0; d<Boundary.size(); d++) this->Cells.push_back(Boundary[d].size());
}
