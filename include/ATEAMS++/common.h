
#ifndef ATEAMS_COMMON_H
#define ATEAMS_COMMON_H

#include <SparseRREF/sparse_mat.h>
#include <phat/boundary_matrix.h>
#include <phat/compute_persistence_pairs.h>
#include <phat/representations/default_representations.h>

#include <set>
#include <vector>
#include <map>

#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <filesystem>

namespace ATEAMS {
	/** Unsigned long int. */
	typedef unsigned long ulong;

	/** Default type used for indexing in SparseRREF matrices. */
	typedef uint32_t index_t;

	/** Default type used for data storage in SparseRREF matrices. */
	typedef ulong data_t;

	/** SparseRREF matrix. */
	typedef SparseRREF::sparse_mat<data_t, index_t> ZpMatrix;

	/** SparseRREF vector. */
	typedef SparseRREF::sparse_vec<data_t, index_t> ZpVector;

	/** SparseRREF field, \f$\mathrm{GF}(p) \cong \Z/p \Z\f$ for \f$p\f$ prime. */
	typedef SparseRREF::field_t Zp;
	
	/** @cond */
	/** SparseRREF index of pivots. */
	typedef SparseRREF::pivot_t<index_t> ZpPivot;

	/** Vector of @ref ATEAMS::ZpMatrix. */
	typedef std::vector<ZpMatrix> ZpMatrices;

	/** Vector of @ref ATEAMS::ZpVector. */
	typedef std::vector<ZpVector> ZpVectors;

	/** Alias of @ref ATEAMS::ZpVectors. */
	typedef std::vector<ZpVector> ZpBasis;

	/** Vector of vector of @ref ATEAMS::ZpBasis. */
	typedef std::vector<ZpBasis> ZpBases;

	/** Flat boundary matrix. */
	typedef std::vector<std::vector<int>> FlatBoundaryMatrix;

	/** PHAT column. */
	typedef std::vector<phat::index> PHATColumn;

	/** PHAT boundary matrix. */
	typedef phat::boundary_matrix<phat::bit_tree_pivot_column> PHATBoundaryMatrix;

	/** PHAT persistence pairs. */
	typedef phat::persistence_pairs PHATPairs;

	/** PHAT persistence algorithm. */
	typedef phat::twist_reduction PHATTwist;
	/** @endcond */


	/**
	 * @class DataWriter
	 * @brief Utility function for writing various data types to file.
	 */
	class DataWriter {
		public:
			/**
			 * @brief Writes a @ref ATEAMS::ZpMatrix to file.
			 * 
			 * @param S Storage.
			 * @param destination Filepath.
			 */
			static inline void write(ATEAMS::ZpMatrix S, std::string destination) {
				checkDestinationExists(destination);

				std::ofstream out;
				out.open(destination);

				SparseRREF::sparse_mat_write<ATEAMS::data_t,std::ofstream,ATEAMS::index_t>(
					S, out, SparseRREF::SPARSE_FILE_TYPE_SMS
				);

				out.close();
			};

			/**
			 * @brief Writes a @ref ATEAMS::ZpMatrix to file.
			 * 
			 * @param S Storage.
			 * @param destination Filepath.
			 */
			static inline void write(std::vector<int> S, std::string destination) {
				checkDestinationExists(destination);

				// Check 
				std::ofstream out;
				out.open(destination);

				for (auto k : S) {
					out << k << std::endl;
				}

				out.close();
			};

			/**
			 * @brief Writes a @ref ATEAMS::ZpMatrix to file.
			 * 
			 * @param S Storage.
			 * @param destination Filepath.
			 */
			static inline void write(std::vector<double> S, std::string destination) {
				checkDestinationExists(destination);

				std::ofstream out;
				out.open(destination);

				for (auto k : S) {
					out << std::format("{:.8f}", k) << std::endl;
				}

				out.close();
			};

			/**
			 * @brief Writes a @ref ATEAMS::ZpMatrix to file.
			 * 
			 * @param S Storage.
			 * @param destination Filepath.
			 */
			static inline void write(std::vector<float> S, std::string destination) {
				checkDestinationExists(destination);

				std::ofstream out;
				out.open(destination);

				for (auto k : S) {
					out << std::format("{:.8f}", k) << std::endl;
				}

				out.close();
			};

		private:
			static inline void checkDestinationExists(std::string destination) {
				std::filesystem::path writable(destination);

				if (!std::filesystem::exists(writable.parent_path())) {
					std::cerr << std::format("Couldn't find filepath {}, forcibly writing...", writable.parent_path().string()) << std::endl;
					std::filesystem::create_directory(writable.parent_path());
				}
			}
	};
}


#endif
