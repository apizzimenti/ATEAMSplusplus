
#ifndef ATEAMS_ARITHMETIC_OPTIONS_H
#define ATEAMS_ARITHMETIC_OPTIONS_H

#include "ATEAMS++/util.h"

#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_rref.h>

#include <omp.h>
#include <thread>


namespace ATEAMS::arithmetic {
	/**
	 * @brief Alias of `SparseRREF::rref_options`.
	 */
	typedef SparseRREF::rref_option RREFOptions;

	/**
	 * @brief Alias of `SparseRREF::rref_options`.
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

	/**
	 * @class ComputeOptions
	 * @brief Convenience class for managing Flint and multithreading. See
	 * 	@ref ATEAMS::statistics::Chain; not something users should need.
	 * 
	 * @var ComputeOptions::opt
	 * 	@brief Pointer to a @ref ATEAMS::arithmetic::ComputeOptions object
	 * 	used to keep track of the thread pool and process-killing keystrokes.
	 * 
	 * @var ComputeOptions::serial
	 * 	@brief Pointer to a @ref ATEAMS::arithmetic::SerialContainers, reusable
	 * 	containers for serial persistence computation.
	 * 
	 * @var ComputeOptions::parallel
	 * 	@brief Pointer to a @ref ATEAMS::arithmetic::ParallelContainers, reusable
	 * 	containers for parallel persistence computation.
	 */
	template <typename RingLike>
	class ComputeOptions {
		public:
			RREFOptions* opt;
			SerialContainers* serial;
			ParallelContainers* parallel;

			/**
			 * @brief Constructor.
			 */
			ComputeOptions() {
				this->opt = new RREFOptionType;

				this->parallel = new ParallelContainers;
				this->parallel->enabled = true;

				this->serial = new SerialContainers;
			};

			/**
			 * @brief Initializes a multithreaded computing environment.
			 * @returns An execution thread.
			 */
			std::thread spinUp() {
				Flint::set_memory_functions();
				this->opt->pool.reset();
				this->opt->method = 0;

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

