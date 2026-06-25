
[![docs](https://img.shields.io/badge/%E2%93%98-Documentation-red)](https://apizzimenti.github.io/ATEAMSplusplus/)
[![arXiv](https://img.shields.io/badge/arXiv-2507.13503-green.svg)](https://arxiv.org/abs/2507.13503)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.20817889-blue)](https://doi.org/10.5281/zenodo.20817889)

# Algebraic Topology-Enabled AlgorithMs for Spin systems in C++ (ATEAMS++)

ATEAMS++ is a high-performance C++ library for simulating generalized Potts, random-cluster, and percolation models at scale. It depends directly on [SparseRREF](https://github.com/munuxi/SparseRREF) and PHAT [[1](https://www.sciencedirect.com/science/article/pii/S0747717116300098),[2](https://bitbucket.org/phat-code/phat/src/master/)] for fast linear algebra over finite fields. It is a C++ port of the original [ATEAMS](https://github.com/apizzimenti/ATEAMS) library.


### Installation

1. **Clone the repository.** Run
	```zsh
	$ git clone --recursive https://github.com/apizzimenti/ATEAMSplusplus.git
	```
	in your favorite shell. The `--recursive` flag ensures that SparseRREF is actually cloned, not just referenced.

2. **Install dependencies.** ATEAMS++ requires

	1. [OpenMP](https://github.com/openmp) and
	2. [FLINT](https://flintlib.org/),

	which are both available through most package managers. [Make sure your compiler complies with the C++23 standard](https://en.cppreference.com/cpp/compiler_support/23); GCC ≥ 15 is a good bet.

3. **Build and install the library.** In your favorite shell, running
	```zsh
	$ cd ATEAMSplusplus
	$ make build
	$ make install
	```
	will build ATEAMS++ to the `build` directory, then install its headers and bundled dependencies (SparseRREF, PHAT) to your system's default install location. If you want to install it somewhere else,
	```zsh
	$ sudo cmake --install build --prefix "somewhere/else"
	```
	will do the trick.

4. **Double-check.** To verify that the library is installed correctly, running
	```zsh
	$ pkg-config --cflags --libs ATEAMS++
	```
	should give something like
	```
	-L/usr/local/lib -lATEAMS++ -lflint -I/usr/local/include -std=c++23
	```

### Usage

Most experiments follow a straightforward template:

1. create a Complex;
2. create a Model that consumes the Complex;
3. run a Chain on the Model.

The below example samples plaquette Bernoulli percolation 1000 times on the scale-6 4-fold torus:

```C++
#include <ATEAMS++/ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

int main() {
	// Construct a cubical complex.
	vector<int> corners = {6,6,6,6};
	complexes::Cubical C(corners, true);

	// Parametrize the model.
	models::BernoulliParameters params;
	params.dimension = 2;

	// Construct the model, then the chain.
	models::Bernoulli percolation(&C, params);

	int iterations = 1000;
	using Chain = statistics::Chain<models::Bernoulli>;
	Chain M(&percolation, iterations);

	// Iterate over the chain, storing the rank of the 2nd homology group as we
	// go.
	vector<int> ranks;

	for (models::BernoulliState* state : M.simulate<models::BernoulliState>()) {
		ranks.push_back(state->rank);
	}

	// Compute the expected rank of the 2nd homology group.
	int T = std::accumulate(ranks.begin(), ranks.end(), 0);
	double e = (double)T/(double)iterations;

	std::cout << "expected rank is " << std::format("{:.2f}", e) << std::endl;
	
	return 0;
}
```

If the above snippet is in a file called `main.cpp`, then compiling the program and running the executable should give something like

```zsh
$ g++ `pkg-config --libs --cflags ATEAMS++` main.cpp -o main
$ ./main
expected rank is 3.09
```

When timed, this takes ~11 seconds.

### Known issues
	* On some Linux systems, linking against ATEAMS++ after it's installed may
	cause some issues. In the event you get a linking error talking about undefined
	symbols in a `tbb` namespace, you may have to recompile your program with
	```zsh
	$ $CXX yourfile.cpp `pkg-config --libs --cflags ATEAMS++` -ltbb -o executable
	```
	for things to work.
