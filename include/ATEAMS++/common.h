
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

	/** Finite-field type for SparseRREF machinery. */
	typedef ulong ff;

	/** Rational type for SparseRREF machinery. */
	typedef SparseRREF::rat_t rational;	

	/** SparseRREF matrix. */
	template <typename T=ff>
	using SparseMatrix = SparseRREF::sparse_mat<T,index_t>;

	/** SparseRREF vector. */
	template <typename T=ff>
	using SparseVector = SparseRREF::sparse_vec<T,index_t>;

	/** SparseRREF field, \f$\mathrm{GF}(p) \cong \Z/p \Z\f$ for \f$p\f$ prime. */
	typedef SparseRREF::field_t Field;
	
	/** @cond */
	/** SparseRREF index of pivots. */
	typedef SparseRREF::pivot_t<index_t> SparsePivot;

	/** Vector of @ref ATEAMS::SparseMatrix. */
	template <typename T=ff>
	using SparseMatrices = std::vector<SparseMatrix<T>>;

	/** Vector of @ref ATEAMS::SparseVectors. */
	template <typename T=ff>
	using SparseVectors = std::vector<SparseVector<T>>;

	/** Alias of @ref ATEAMS::SparseVectors. */
	template <typename T=ff>
	using SparseBasis = std::vector<SparseVector<T>>;

	/** Vector of vector of @ref ATEAMS::SparseBasis. */
	template <typename T>
	using SparseBases = std::vector<SparseBasis<T>>;

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


	// /**
	//  * @class DataWriter
	//  * @brief Utility function for writing various data types to file.
	//  */
	// class DataWriter {
	// 	public:
	// 		/**
	// 		 * @brief Writes a @ref ATEAMS::ZpMatrix to file.
	// 		 * 
	// 		 * @param S Storage.
	// 		 * @param destination Filepath.
	// 		 */
	// 		static inline void write(ATEAMS::ZpMatrix S, std::string destination) {
	// 			checkDestinationExists(destination);

	// 			std::ofstream out;
	// 			out.open(destination);

	// 			SparseRREF::sparse_mat_write<ATEAMS::data_t,std::ofstream,ATEAMS::index_t>(
	// 				S, out, SparseRREF::SPARSE_FILE_TYPE_SMS
	// 			);

	// 			out.close();
	// 		};

	// 		/**
	// 		 * @brief Writes a @ref ATEAMS::ZpMatrix to file.
	// 		 * 
	// 		 * @param S Storage.
	// 		 * @param destination Filepath.
	// 		 */
	// 		static inline void write(std::vector<int> S, std::string destination) {
	// 			checkDestinationExists(destination);

	// 			// Check 
	// 			std::ofstream out;
	// 			out.open(destination);

	// 			for (auto k : S) {
	// 				out << k << std::endl;
	// 			}

	// 			out.close();
	// 		};

	// 		/**
	// 		 * @brief Writes a @ref ATEAMS::ZpMatrix to file.
	// 		 * 
	// 		 * @param S Storage.
	// 		 * @param destination Filepath.
	// 		 */
	// 		static inline void write(std::vector<double> S, std::string destination) {
	// 			checkDestinationExists(destination);

	// 			std::ofstream out;
	// 			out.open(destination);

	// 			for (auto k : S) {
	// 				out << std::format("{:.8f}", k) << std::endl;
	// 			}

	// 			out.close();
	// 		};

	// 		/**
	// 		 * @brief Writes a @ref ATEAMS::ZpMatrix to file.
	// 		 * 
	// 		 * @param S Storage.
	// 		 * @param destination Filepath.
	// 		 */
	// 		static inline void write(std::vector<float> S, std::string destination) {
	// 			checkDestinationExists(destination);

	// 			std::ofstream out;
	// 			out.open(destination);

	// 			for (auto k : S) {
	// 				out << std::format("{:.8f}", k) << std::endl;
	// 			}

	// 			out.close();
	// 		};

	// 	private:
	// 		static inline void checkDestinationExists(std::string destination) {
	// 			std::filesystem::path writable(destination);

	// 			if (!std::filesystem::exists(writable.parent_path())) {
	// 				std::cerr << std::format("Couldn't find filepath {}, forcibly writing...", writable.parent_path().string()) << std::endl;
	// 				std::filesystem::create_directory(writable.parent_path());
	// 			}
	// 		}
	// };
}


#endif
