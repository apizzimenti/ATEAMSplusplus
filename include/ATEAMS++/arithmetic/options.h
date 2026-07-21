
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
	class ComputeOptions {
		public:
			RREFOptions* opt;

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
	};
}


#endif

