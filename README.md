
[![arXiv](https://img.shields.io/badge/arXiv-2507.13503-green.svg)](https://arxiv.org/abs/2507.13503)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.20817889-blue)](https://doi.org/10.5281/zenodo.20817889)
[![docs](https://img.shields.io/badge/%E2%93%98-Documentation-red)](https://apizzimenti.github.io/ATEAMSplusplus/)
[![License:MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Algebraic Topology-Enabled AlgorithMs for Spin systems, in C++

ATEAMS++ is a high-performance C++ library for simulating
[generalized Potts, random-cluster, and percolation models](https://arxiv.org/abs/2507.13503) at scale. It depends directly on [SparseRREF](https://github.com/munuxi/SparseRREF) and PHAT [[1](https://www.sciencedirect.com/science/article/pii/S0747717116300098),[2](https://bitbucket.org/phat-code/phat/src/master/)] for fast linear algebra over finite fields and [`BS::thread_pool`](https://github.com/bshoshany/thread-pool) for high-control parallelism. It is a C++ port of the original [ATEAMS](https://github.com/apizzimenti/ATEAMS) library, now deprecated in favor of this software.


## Prerequisites

> [!TIP]
> The [Homebrew](https://brew.sh/) package manager is strongly recommended for users running macOS, any flavor of Linux, or Windows Subsystem for Linux (WSL).

The only real prerequisite is **a compiler that supports the C++23 standard.** GCC 15 is the oldest cross-compatible compiler that supports all ATEAMS++'s features. Newer versions of GCC (e.g. 16) work as well.

* **on macOS**,
	Install [GCC 15](https://formulae.brew.sh/formula/gcc@15) [(or 16)](https://formulae.brew.sh/formula/gcc) through Homebrew. After installation, test that you have the correct version by running
	```bash
	$ g++-15 --version

	g++-15 (Homebrew GCC 15.3.0) 15.3.0
	Copyright (C) 2025 Free Software Foundation, Inc.
	This is free software; see the source for copying conditions.  There is NO
	warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	
	```
	You can set `gcc-15` and `g++-15` to be your user's default C/C++ compilers by running
	```bash
	$ echo "export CC=\"gcc-15\" \nexport CXX=\"g++-15\" \nalias gcc=\$CC \nalias g++=\$CXX" >> ~/.zshrc && source ~/.zshrc
	```
	which appends aliases to your zsh configuration file (and creates the file if it doesn't already exist, which it should). You can verify that the aliases were correctly installed by running
	```bash
	$ cat ~/.zshrc

	export CC="gcc-15"
	export CXX="g++-15"
	alias gcc=$CC
	alias g++=$CXX
	```

	You can then run `g++` (with no arguments) just to check:

	```bash
	$ g++
	g++-15: fatal error: no input files
	compilation terminated.
	```

* **on Linux**,
	(Ubuntu 22.04/24.04), use the same [installation steps as on macOS](#macOS) (adjusting your shell if you're using bash), or by following [these instructions](https://askubuntu.com/a/1518433) (adjusting the version from 14 to 15).

* **on Windows**,
	[Try MinGW](https://www.mingw-w64.org/getting-started/msys2/); it should come with a recent version of GCC.
		


## Installation

1. **Clone the repository.** In your favorite shell, navigate to a sensible directory and run
	```bash
	$ git clone --recursive https://github.com/apizzimenti/ATEAMSplusplus.git
	```
	The `--recursive` flag ensures that SparseRREF is actually cloned, not just referenced.

2. **Install dependencies.** ATEAMS++ requires

	1. [OpenMP](https://github.com/openmp) (which comes prepackaged with GCC compilers),
	2. [FLINT](https://flintlib.org/),
	3. [TBB](https://github.com/uxlfoundation/onetbb), and
	4. (optionally) [jemalloc](https://jemalloc.net/),

	which are all available through most package managers. [Make sure your compiler complies with the C++23 standard](https://en.cppreference.com/cpp/compiler_support/23); GCC ≥ 15 is a good bet.

	To install the dependencies...
	* **on macOS**,
		```bash
		$ brew install flint tbb jemalloc
		```
	* **on Linux (Ubuntu)**,
		```bash
		$ sudo apt update
		$ sudo apt install libgmp-dev libmpfr-dev libflint-dev libtbb-dev libjemalloc-dev
		```
	* **on Windows** (using [vcpkg](https://vcpkg.io/en/)),
		```bash
		$ vcpkg install flint tbb jemalloc
		```


> [!NOTE]
> **jemalloc is _not_ required to use this software, but it is recommended:** after extensive testing on machines with 4 (macOS M2), 36 (Intel Xeon W-2295), and 192 (AMD Ryzen Threadripper PRO 7995WX) logical cores, **jemalloc provides the best multi-threaded performance** when compared to no alternate memory allocator, mimalloc, and tcmalloc. If you don't want to use it, set the `jemalloc` variable in the CMake config file (`CMakeLists.txt`) to `FALSE` and it will be excluded.
>
> If you want to use another memory allocator like [mimalloc](https://github.com/microsoft/mimalloc) or [tcmalloc](https://github.com/google/tcmalloc), install it in the same way as jemalloc:
> * **on macOS**,
> 	```bash
> 	$ brew install gperftools mimalloc
> 	```
> * **on Linux (Ubuntu)**,
> 	```bash
> 	$ sudo apt update
> 	$ sudo apt install libgoogle-perftools-dev libmimalloc-dev
> 	```
> * **on Windows** (using [vcpkg](https://vcpkg.io/en/)),
> 	```bash
> 	$ vcpkg install gperftools mimalloc
>	```
> You can set any of the `jemalloc`, `tcmalloc`, or `mimalloc` tools to `TRUE` in the CMake config file (`CMakeLists.txt`), but the build will fail if more than one of the variables is set to `TRUE`.

</br>

3. **Build and install the library.** In your favorite shell, running
	```bash
	$ cd ATEAMSplusplus
	$ make build test install
	```
	will build ATEAMS++'s tests to the `build` directory, run those tests, then install its headers and bundled dependencies (SparseRREF, PHAT) to your system's default install location. If you want to install it somewhere else,
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
	-I/usr/local/include -std=c++23 -fopenmp -Ofast -march=native -flto -funroll-loops -fprefetch-loop-arrays -fomit-frame-pointer -fno-rtti -L/usr/local/lib -lflint -L/opt/homebrew/lib -ltbb -ljemalloc
	```

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

// Function for timing the loop; returns fractional seconds.
double timer(auto loop) {
	auto start = chrono::steady_clock::now();
	loop();
	auto end = chrono::steady_clock::now();

	auto delta = chrono::duration_cast<chrono::milliseconds>(end-start);
	return delta.count()/1000.0;
}

int main() {
	int iterations = 1000;
	int dimension = 2;

	// Construct a 6x6x6x6 cubical 4-torus, initialize a 2-dimensional Bernoulli
	// percolation model on the torus, and initialize a Markov chain to simulate
	// the percolation model.
	complexes::Cubical<Z2> plex({6,6,6,6});
	models::Bernoulli percolation(&plex, dimension);
	statistics::Chain<models::Bernoulli> M(&percolation, iterations);

	// Create a vector to store the ranks of the (dimension)th homology groups
	// of each subcomplex sampled by the model.
	vector<int> ranks(M.steps, 0);

	// Time the loop.
	double elapsed = timer([&]{
		for (models::Bernoulli::State state : M.simulate()) {
			ranks[state.t] = state.rank;
		}
	});

	// Compute the expected rank.
	double e = accumulate(ranks.begin(), ranks.end(), 0)/(double)M.steps;

	// Pretty-print.
	cout << std::format(
		"expected rank is {:.3f} :: took {:.3f} seconds :: ATEAMS++ v{}",
		e,
		elapsed,
		ATEAMS_VERSION
	) << endl;

	return 0;
}
```

If the above snippet is in a file called `main.cpp`, then compiling the program and running the executable should give something like

```bash
$ $CXX main.cpp `pkg-config --libs --cflags ATEAMS++` -o main
$ ./main
expected rank is 3.09
```

This takes ~8 seconds on my 2022 MacBook Air with an Apple M2 processor.

> [!TIP]
> If you plan to run large-scale simulations on an external machine, the [`magnetizationplusplus` workflow](https://github.com/apizzimenti/magnetizationplusplus) may be of use to you.


## Known issues
* `InvadedCluster` and `Invasion` perform $`\mathbb Z/2\mathbb Z`$ persistence when using rational coefficients.
* Some models don't have complete debugging information.


## Contributing
We welcome contributions. A few housekeeping rules:

1. **Use the fork + branch + PR flow.** Fork + clone this repository, then create a new branch for each new feature you add. Once it's camera-ready, open a pull request on the main repository page, and we'll review it.
2. **Test!** Please write unit tests (like the ones in `tests/`) for any significant features (e.g. new Models, arithmetic functions, parallelization) you add. **_Pull requests without accompanying unit tests will be automatically closed._**

## Citing

```bibtex
@software{ATEAMS++,
	title={{ATEAMS++: Algebraic Topology-Enabled AlgorithMs for Spin systems, in C++}},
	author={Pizzimenti, Anthony E. and Eldridge, Summer and Schweinhart, Ben and Duncan, Paul},
	url={github.com/apizzimenti/ATEAMSplusplus},
	version={1.0.0},
	doi={10.5281/zenodo.20817889}
}
```


