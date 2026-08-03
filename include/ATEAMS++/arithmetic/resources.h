
#ifndef ATEAMS_ARITHMETIC_OPTIONS_H
#define ATEAMS_ARITHMETIC_OPTIONS_H

#include "ATEAMS++/util.h"

#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_rref.h>

#include <omp.h>
#include <thread>


namespace ATEAMS::arithmetic {
	/**
	 * @brief Alias of `SparseRREF::rref_resources`.
	 */
	typedef SparseRREF::rref_option RREFOptions;

	/**
	 * @brief Alias of `SparseRREF::rref_resources`.
	 */
	typedef SparseRREF::rref_option_t RREFOptionType;

	/**
	 * @struct SerialContainers
	 * @brief Reusable containers for serial persistence computation. See
	 * @ref ATEAMS::topology::persistence::twist.
	 * 
	 * @var SerialContainers::marked
	 * 	Set for marking cycle-inducing cells.
	 * 
	 * @var SerialContainers::lookup
	 * 	Lookup table for youngest chains sharing a face.
	 */
	struct SerialContainers {
		public:
			std::set<int> marked;
			std::vector<int> lookup;

		/**
		 * @brief Constructs the lookup table.
		 * 
		 * @param cells Number of cells in the complex.
		 */
		void build(int cells) {
			this->lookup = std::vector<int>(cells, 0);
		}

		/**
		 * @brief Flushes the containers for reuse.
		 */
		void flush() {
			this->marked.clear();
			std::fill(this->lookup.begin(), this->lookup.end(), 0);
			// for (int t=0; t < this->lookup.size(); t++) this->lookup[t] = 0;
		}
	};

	/**
	 * @struct ParallelContainers
	 * @brief Reusable containers for parallel persistence computation. See
	 * @ref ATEAMS::topology::persistence::JIT.
	 * 
	 * @var ParallelContainers::marked
	 * 	Vector of sets for marking cycle-inducing cells. Created as a vector to
	 * 	ensure all threads write to different locations.
	 * 
	 * @var ParallelContainers::zeroed
	 * 	Allows (e.g.) @ref ATEAMS::topology::persistence::JIT to simulate the
	 * 	clearing optimization of the twist algorithm.
	 * 
	 * @var ParallelContainers::lookup
	 * 	Lookup table for youngest chains sharing a face.
	 * 
	 * @var ParallelContainers::enabled
	 * 	Is parallel computing enabled? Defaults to `true`.
	 */
	struct ParallelContainers {
		public:
			std::vector<std::set<int>> marked;
			std::vector<int> lookup;
			std::vector<bool> zeroed;
			bool enabled = true;

			/**
			 * @brief Constructs data structures.
			 * 
			 * @param cells Number of cells in the complex.
			 * @param blocks Number of matrix blocks (dimensions) in the complex.
			 */
			void build(int cells, int blocks) {
				this->zeroed = std::vector<bool>(cells, false);
				this->lookup = std::vector<int>(cells, 0);
				this->marked = std::vector<std::set<int>>(blocks, std::set<int>());
			}

			/**
			 * @brief Flushes data structures for reuse.
			 */
			void flush() {
				for (int t=0; t < this->marked.size(); t++) this->marked[t].clear();
				std::fill(this->lookup.begin(), this->lookup.end(), 0);
				std::fill(this->zeroed.begin(), this->zeroed.end(), false);
			}
	};

	template <typename RingLike>
	struct FiniteArithmetic {
		std::vector<std::vector<typename RingLike::dtype>> add;
		std::vector<std::vector<typename RingLike::dtype>> multiply;
		std::vector<typename RingLike::dtype> negate;
		std::vector<typename RingLike::dtype> invert;

		void arithmetize(Ring* R) {
			// Create addition table.
			int c = R->characteristic;
			std::vector<std::vector<typename RingLike::dtype>> addition(c, std::vector<typename RingLike::dtype>(c, 0));

			for (int a=0; a < c; a++) {
				for (int b=0; b < c; b++) {
					addition[a][b] = (typename RingLike::dtype)((a+b) % c);
				}
			}

			// Create multiplication table.
			std::vector<std::vector<typename RingLike::dtype>> multiplication(c, std::vector<typename RingLike::dtype>(c, 0));
			for (int a=0; a < c; a++) {
				for (int b=0; b < c; b++) {
					multiplication[a][b] = (typename RingLike::dtype)((a*b) % c);
				}
			}

			// Create negation table.
			std::vector<typename RingLike::dtype> negation(c, 0);
			for (int a=0; a < c; a++) {
				negation[a] = (typename RingLike::dtype)((c-a)%c);
			}

			// Create inversion table.
			std::vector<typename RingLike::dtype> inversion(c, 0);
			for (int a=0; a < c; a++) {
				inversion[a] = (typename RingLike::dtype)(((int)pow(a, c-2))%c);
			}

			this->add = addition;
			this->multiply = multiplication;
			this->negate = negation;
			this->invert = inversion;
		};
	};

	/**
	 * @class ComputeResources
	 * @brief Convenience class for managing Flint and multithreading. See
	 * 	@ref ATEAMS::statistics::Chain; not something users should need.
	 * 
	 * @var ComputeResources::opt
	 * 	@brief Pointer to a @ref ATEAMS::arithmetic::ComputeResources object
	 * 	used to keep track of the thread pool and process-killing keystrokes.
	 * 
	 * @var ComputeResources::serial
	 * 	@brief Pointer to a @ref ATEAMS::arithmetic::SerialContainers, reusable
	 * 	containers for serial persistence computation.
	 * 
	 * @var ComputeResources::parallel
	 * 	@brief Pointer to a @ref ATEAMS::arithmetic::ParallelContainers, reusable
	 * 	containers for parallel persistence computation.
	 */
	template <typename RingLike>
	class ComputeResources {
		public:
			RREFOptions* opt;
			SerialContainers* serial;
			ParallelContainers* parallel;
			FiniteArithmetic<RingLike>* arithmetic;


			/**
			 * @brief Constructor.
			 */
			ComputeResources() {
				this->opt = new RREFOptionType;

				this->parallel = new ParallelContainers;
				this->parallel->enabled = true;

				this->serial = new SerialContainers;

				this->arithmetic = new FiniteArithmetic<RingLike>;
			};

			void arithmetize(Ring* R) {
				this->arithmetic->arithmetize(R);
			}

			/**
			 * @brief Initializes a multithreaded computing environment.
			 * @returns An execution thread.
			 */
			std::thread spinUp() {
				Flint::set_memory_functions();
				this->opt->pool.reset();
				this->opt->method = 0;
				this->opt->verbose = true;

				return std::thread(key_listener, std::ref(this->opt->abort));
			};

			/**
			 * @brief "Spins down" the multithreaded computing environment.
			 * 
			 * @param listener Pointer to a `std::thread`.
			 */
			void spinDown(std::thread* listener) {
				this->opt->abort = true;
				Flint::clear_cache();
				listener->join();
			};
	};
}


#endif

