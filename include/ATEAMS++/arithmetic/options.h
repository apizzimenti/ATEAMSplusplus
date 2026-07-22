
#ifndef ATEAMS_ARITHMETIC_OPTIONS_H
#define ATEAMS_ARITHMETIC_OPTIONS_H

#include "ATEAMS++/util.h"

#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_rref.h>

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

	/**
	 * @class ParallelOptions
	 * @brief Facilities for more efficient parallel computing.
	 */
	template <typename RingLike>
	struct ParallelOptions {
		public:
			int threads;
			bool enabled = false;
			std::vector<std::vector<int>> indexBlocks;
			std::vector<SparseVector<RingLike>> lScratch;
			std::vector<SparseVector<RingLike>> rScratch;
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
			ParallelOptions<RingLike>* parallel;

			/**
			 * @brief Constructor.
			 */
			ComputeOptions() {
				this->opt = new RREFOptionType;
			};

			/**
			 * @brief Initializes a multithreaded computing environment.
			 * @returns An execution thread.
			 */
			std::thread spinUp() {
				Flint::set_memory_functions();
				this->opt->pool.reset();
				this->opt->method = 0;

				this->initializeParallelism();

				return std::thread(key_listener, std::ref(this->opt->abort));
			};

			/**
			 * @brief Initializes a multithreaded computing environment.
			 * @returns An execution thread.
			 */
			std::thread spinUp(int strategy) {
				Flint::set_memory_functions();
				this->opt->pool.reset();
				this->opt->method = 0;

				this->initializeParallelism(strategy);

				return std::thread(key_listener, std::ref(this->opt->abort));
			};

			/**
			 * @brief "Spins down" the multithreaded computing environment.
			 */
			void spinDown(std::thread* listener) {
				this->opt->abort = true;
				Flint::clear_cache();
				listener->join();
				// listener->detach();
			};

		private:
			void initializeParallelism(int strategy) {
				// Determine how many threads we're using; this probably won't
				// change during execution.
				this->parallel = new ParallelOptions<RingLike>;
				this->parallel->enabled = true;

				int threads = this->opt->pool.get_thread_count();
				this->parallel->threads = std::max(1, (int)(pow(2, -strategy)*threads));
				this->opt->pool.reset(this->parallel->threads);

				this->parallel->indexBlocks = std::vector<std::vector<int>>(this->parallel->threads, std::vector<int>(2,0));
				this->parallel->lScratch = std::vector<SparseVector<RingLike>>(this->parallel->threads);
				this->parallel->rScratch = std::vector<SparseVector<RingLike>>(this->parallel->threads);
			};

			void initializeParallelism() {
				// Determine how many threads we're using; this probably won't
				// change during execution.
				this->parallel = new ParallelOptions<RingLike>;
				this->parallel->enabled = true;
				
				int threads = this->opt->pool.get_thread_count();
				this->parallel->threads = 2*threads+1;

				this->parallel->indexBlocks = std::vector<std::vector<int>>(this->parallel->threads, std::vector<int>(2,0));
				this->parallel->lScratch = std::vector<SparseVector<RingLike>>(this->parallel->threads);
				this->parallel->rScratch = std::vector<SparseVector<RingLike>>(this->parallel->threads);
			};
	};
}


#endif

