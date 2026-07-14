
[![arXiv](https://img.shields.io/badge/arXiv-2507.13503-green.svg)](https://arxiv.org/abs/2507.13503)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.20817889-blue)](https://doi.org/10.5281/zenodo.20817889)
[![docs](https://img.shields.io/badge/%E2%93%98-Documentation-red)](https://apizzimenti.github.io/ATEAMSplusplus/)
[![License:MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Algebraic Topology-Enabled AlgorithMs for Spin systems, in C++

ATEAMS++ is a high-performance C++ library for simulating
[generalized Potts, random-cluster, and percolation models](https://arxiv.org/abs/2507.13503) at scale. It depends directly on [SparseRREF](https://github.com/munuxi/SparseRREF) and PHAT [[1](https://www.sciencedirect.com/science/article/pii/S0747717116300098),[2](https://bitbucket.org/phat-code/phat/src/master/)] for fast linear algebra over finite fields. It is a C++ port of the original [ATEAMS](https://github.com/apizzimenti/ATEAMS) library, which is now deprecated in favor of this software.


## Installation

1. **Clone the repository.** Run
	```bash
	$ git clone --recursive https://github.com/apizzimenti/ATEAMSplusplus.git
	```
	in your favorite shell. The `--recursive` flag ensures that SparseRREF is actually cloned, not just referenced.

2. **Install dependencies.** ATEAMS++ requires

	1. [OpenMP](https://github.com/openmp),
	2. [FLINT](https://flintlib.org/), and
	3. (optionally!) [mimalloc](https://github.com/microsoft/mimalloc),

	which are all available through most package managers. [Make sure your compiler complies with the C++23 standard](https://en.cppreference.com/cpp/compiler_support/23); GCC ≥ 15 is a good bet. mimalloc is _not_ required to use this software, but it is recommended; if you don't want to use it, set the `mimalloc` variable in `CMakeLists.txt` to `FALSE` and it will be excluded.

3. **Build and install the library.** In your favorite shell, running
	```bash
	$ cd ATEAMSplusplus
	$ make build
	$ make install
	```
	will build ATEAMS++ to the `build` directory, then install its headers and bundled dependencies (SparseRREF, PHAT) to your system's default install location. If you want to install it somewhere else,
	```bash
	$ sudo cmake --install build --prefix "somewhere/else"
	```
	will do the trick.

4. **Double-check.** To verify that the library is installed correctly, running
	```bash
	$ pkg-config --cflags --libs ATEAMS++
	```
	should give something like
	```bash
	-I/usr/local/include -std=c++23 -O3 -L/usr/local/lib -lflint
	```

5. **Triple-check.** If you want to run pre-built tests, set the `TEST` variable in `CMakeLists.txt` to `TRUE`, á la
	```
	set(TEST TRUE)
	```
	then run `make test`.

## Usage

Most experiments follow a straightforward template:

1. create a Complex;
2. create a Model that consumes the Complex;
3. run a Chain on the Model.

The below example samples $`2`$-dimensional plaquette Bernoulli percolation 1,000 times on the scale-$`6`$ $`4`$-fold torus $`\mathbb T_6^4`$:

```cpp
#include <ATEAMS++.h>

using namespace ATEAMS;
using namespace std;

int main() {
	int iterations = 1000;
	int dimension = 2;

	complexes::Cubical<Z2> plex({6,6,6,6});
	models::Bernoulli percolation(&plex, dimension);
	statistics::Chain<models::Bernoulli> M(&percolation, iterations);

	vector<int> ranks(M.steps, 0);

	for (models::Bernoulli::State state : M.simulate()) { ranks[state.t] = state.rank; }

	double e = (double)accumulate(ranks.begin(), ranks.end(), 0)/(double)M.steps;
	cout << "expected rank is " << format("{:.2f}", e) << endl;

	return 0;
}
```

If the above snippet is in a file called `main.cpp`, then compiling the program and running the executable should give something like

```bash
$ $CXX main.cpp `pkg-config --libs --cflags ATEAMS++` -o main
$ ./main
expected rank is 3.09
```

This takes ~8 seconds on a 2022 MacBook Air with an Apple M2 processor. **If you plan to run large-scale simulations on an external machine, the [`magnetizationplusplus` workflow](https://github.com/apizzimenti/magnetizationplusplus) may be of use to you.**

## Known issues
* On some Linux systems, linking against ATEAMS++ after it's installed may cause some issues. In the event you get a linking error talking about undefined symbols in a `tbb` namespace, you may have to recompile your program with

	```bash
	$ $CXX yourfile.cpp `pkg-config --libs --cflags ATEAMS++` -ltbb -o executable
	```

	for things to work.

* Using `mimalloc` on macOS with GCC ≥ 16 can sometimes yield slower compute times and segfaults.
* `InvadedCluster` and `Invasion` perform $`\mathbb Z/2\mathbb Z`$ persistence when using rational coefficients.
* Some models don't have complete debugging information.


## Contributing
We welcome contributions. A few housekeeping rules:

1. **Use the fork + branch + PR flow.** Fork + clone this repository, then create a new branch for each new feature you add. Once it's camera-ready, open a pull request on the main repository page, and we'll review it.
2. **Test!** Please write unit tests (like the ones in `tests/`) for any significant features (e.g. new Models, arithmetic functions, parallelization) you add. _Pull requests without accompanying unit tests will be automatically closed._

## Citing

```bibtex
@software{ATEAMS++,
	title={{ATEAMS++: Algebraic Topology-Enabled AlgorithMs for Spin systems, in C++}},
	author={Pizzimenti, Anthony E. and Eldridge, Summer and Schweinhart, Ben and Duncan, Paul},
	url={github.com/apizzimenti/ATEAMSplusplus},
	version={1.8},
	doi={10.5281/zenodo.20817889}
}
```


