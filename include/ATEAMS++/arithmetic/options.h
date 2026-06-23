
#ifndef ATEAMS_ARITHMETIC_OPTIONS_H
#define ATEAMS_ARITHMETIC_OPTIONS_H

#include "ATEAMS++/util.h"
#include <SparseRREF/sparse_mat.h>
#include <SparseRREF/sparse_rref.h>

#include <thread>


namespace ATEAMS {
	namespace arithmetic {
		/**
		 * @brief Convenience renaming of `SparseRREF::rref_options`.
		 */
		typedef SparseRREF::rref_option ComputeOptions;

		/**
		 * @brief Convenience renaming of `SparseRREF::rref_options`.
		 */
		typedef SparseRREF::rref_option_t ComputeOption;

		/**
		 * @brief Convenience class for managing Flint and multithreading.
		 */
		class ThreadOptions {
			public:
				/**
				 * @brief (Pointer to) a `SparseRREF::rref_option` object; used
				 * 	to throw flags when computations are complete.
				 */
				ComputeOptions* opt;

				ThreadOptions() {
					this->opt = new ComputeOption;
				};

				/**
				 * @brief "Spins up" the appropriate multithreaded computing environment.
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
}


#endif

