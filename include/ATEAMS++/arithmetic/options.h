
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
	typedef SparseRREF::rref_option ComputeOptions;

	/**
	 * @brief Alias of `SparseRREF::rref_options`.
	 */
	typedef SparseRREF::rref_option_t ComputeOption;

	/**
	 * @class ThreadOptions
	 * @brief Convenience class for managing Flint and multithreading. See
	 * 	@ref ATEAMS::statistics::Chain; not something users should need.
	 * 
	 * @var ThreadOptions::opt
	 * 	@brief Pointer to a @ref ATEAMS::arithmetic::ComputeOptions object
	 * 	used to keep track of the thread pool and process-killing keystrokes.
	 * 
	 * @var ThreadOptions::parallelSparseAddition
	 * 	@brief Enables parallelization for @ref ATEAMS::arithmetic::SparseVectorAddition.
	 * 	Default `true`.
	 */
	class ThreadOptions {
		public:
			ComputeOptions* opt;
			bool parallelSparseAddition = true;

			/**
			 * @brief Constructor.
			 */
			ThreadOptions() {
				this->opt = new ComputeOption;
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

