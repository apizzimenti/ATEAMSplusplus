
#include "ATEAMS++/complexes/Cubical.h"

#include <cmath>
#include <algorithm>

using namespace std;

/**
 * Gets lexigraphic-order combinations of elements of a vector.
 */
template <typename T>
void comb(vector<T> in, vector<vector<T>> &combs, vector<T> &active, int r, int cursor, int depth) {
	for (int i=cursor; i<in.size()-(r-depth-1); i++) {
		active[depth] = in[i];
		if (depth < r-1) comb(in, combs, active, r, i+1, depth+1);
		else combs.push_back(active);
	}
}

template <typename T>
vector<vector<T>> combinations(vector<int> in, int r) {
	vector<vector<T>> c;
	vector<T> active(r);
	comb(in, c, active, r, 0, 0);

	return c;
}


/**
 * Creates (the boundary matrices of) a D-dimensional Hamming cube.
 */

typedef uint16_t bitstring;
typedef vector<vector<vector<bitstring>>> HammingCube;

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


typedef vector<vector<vector<int>>> HammingCubeBoundary;

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

Cubical::Cubical(vector<int> corners, bool periodic=true) {
	this->corners = corners;
	this->periodic = periodic;

	// Create a Hamming cube of the appropriate dimension.
	HammingCube cube = hamming(corners.size());
	HammingCubeBoundary boundary = hammingBoundary(cube);

	for (auto dimension : boundary) {
		for (auto cells : dimension) {
			printvector(cells);
		}
		cout << endl;
	}
}
