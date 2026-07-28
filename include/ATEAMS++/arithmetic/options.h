
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

	// enum class Strategy { MIXED, PARALLEL, SERIAL };

	// /**
	//  * @class ParallelOptions
	//  * @brief Facilities for more efficient parallel computing.
	//  */
	// template <typename RingLike>
	// struct ParallelOptions {
	// 	public:
	// 		int threads;
	// 		bool enabled = false;

	// 		std::vector<std::set<int>> marked;

	// 		void build(int blocks, int length) {
	// 			this->marked = std::vector<std::set<int>>(blocks, std::set<int>());
	// 		}

	// 		void flush() {
	// 			for (int b=0; b < this->marked.size(); b++) this->marked[b].clear();
	// 		};
	// };

	struct SerialContainers {
		public:
			std::set<int> marked;
			std::vector<int> lookup;

		void build(int cells) {
			this->lookup = std::vector<int>(cells, 0);
		}

		void flush() {
			this->marked.clear();
			std::fill(this->lookup.begin(), this->lookup.end(), 0);
			// for (int t=0; t < this->lookup.size(); t++) this->lookup[t] = 0;
		}
	};

	struct ParallelContainers {
		public:
			std::vector<std::set<int>> marked;
			std::vector<int> lookup;
			std::vector<bool> zeroed;
			bool enabled = true;

			void build(int cells, int blocks) {
				this->zeroed = std::vector<bool>(cells, false);
				this->lookup = std::vector<int>(cells, 0);
				this->marked = std::vector<std::set<int>>(blocks, std::set<int>());
			}

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
	 * @var ComputeOptions::parallelSparseAddition
	 * 	@brief Enables parallelization for @ref ATEAMS::arithmetic::SparseVectorAddition.
	 * 	Default `true`.
	 * 
	 * @var ComputeOptions::parallelSparseAdditionChunkWidth
	 * 	@brief The width a chunk needs to be before using multiple threads. Default
	 * 	`512`.
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
			 */
			void spinDown(std::thread* listener) {
				this->opt->abort = true;
				Flint::clear_cache();
				listener->join();
			};
	};
}


#endif

