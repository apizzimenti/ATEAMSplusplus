
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
	typedef ulong FINITE;

	/** Rational type for SparseRREF machinery. */
	typedef SparseRREF::rat_t RATIONAL;

	typedef SparseRREF::field_t Field;

	struct Ring {
		public:
			Field ring;
	};

	/**
	 * @brief Rationals \f$\Q\f$.
	 * 
	 * @var Q::ring SparseRREF ring underlying this Ring.
	 */
	struct Q : public Ring {
		/**
		 * @brief Data type for this ring.
		 */
		typedef RATIONAL dtype;

		/**
		 * @brief Constructor.
		 * 
		 * @param characteristic Characteristic of the field; superfluous.
		 */
		Q(int characteristic) {
			this->ring = Field(SparseRREF::FIELD_QQ);
		};

		/** 
		 * @brief Constructor.
		*/
		Q() {
			this->ring = Field(SparseRREF::FIELD_QQ);
		};
	};

	struct Zp : public Ring {
		typedef FINITE dtype;

		Zp(int characteristic) {
			this->ring = Field(SparseRREF::FIELD_Fp, characteristic);
		};
	};

	// /** SparseRREF matrix. */
	// template <typename T>
	// using SparseMatrix = SparseRREF::sparse_mat<T,index_t>;

	// /** SparseRREF vector. */
	// template <typename T>
	// using SparseVector = SparseRREF::sparse_vec<T,index_t>;

	/** Generic  */
	template <typename RingLike>
	using SparseMatrix = SparseRREF::sparse_mat<typename RingLike::dtype,index_t>;

	template <typename RingLike>
	using SparseVector = SparseRREF::sparse_vec<typename RingLike::dtype,index_t>;

	/** SparseRREF field: either \f$\Q\f$ or \f$\GF(p) \cong \Z/p \Z\f$ for \f$p\f$ prime. */
	typedef SparseRREF::field_t Field;
	// typedef SparseRREF::field_t Ring;

	/** Template alias for `std::vector`. */
	template <typename T>
	using DenseVector = std::vector<T,std::allocator<T>>;

	/** Vector of @ref ATEAMS::SparseMatrix. */
	template <typename RingLike>
	using SparseMatrices = std::vector<SparseMatrix<RingLike>>;

	/** Vector of @ref ATEAMS::SparseVectors. */
	template <typename RingLike>
	using SparseVectors = std::vector<SparseVector<RingLike>>;

	/** Alias of @ref ATEAMS::SparseVectors. */
	template <typename RingLike>
	using SparseBasis = std::vector<SparseVector<RingLike>>;

	/** Vector of vector of @ref ATEAMS::SparseBasis. */
	template <typename RingLike>
	using SparseBases = std::vector<SparseBasis<RingLike>>;

	/** @cond */
	/** Flat boundary matrix. */
	typedef std::vector<std::vector<int>> FlatBoundaryMatrix;

	/** SparseRREF index of pivots. */
	typedef SparseRREF::pivot_t<index_t> SparsePivot;

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
			 * @brief Writes a @ref ATEAMS::SparseMatrix to file.
			 * 
			 * @param S Storage.
			 * @param destination Filepath.
			 */
			static inline void write(ATEAMS::SparseMatrix<Q> S, std::string destination) {
				checkDestinationExists(destination);

				std::ofstream out;
				out.open(destination);

				SparseRREF::sparse_mat_write<Q::dtype,std::ofstream,ATEAMS::index_t>(
					S, out, SparseRREF::SPARSE_FILE_TYPE_SMS
				);

				out.close();
			};

			/**
			 * @brief Writes a @ref ATEAMS::SparseMatrix to file.
			 * 
			 * @param S Storage.
			 * @param destination Filepath.
			 */
			static inline void write(ATEAMS::SparseMatrix<Zp> S, std::string destination) {
				checkDestinationExists(destination);

				std::ofstream out;
				out.open(destination);

				SparseRREF::sparse_mat_write<Zp::dtype,std::ofstream,ATEAMS::index_t>(
					S, out, SparseRREF::SPARSE_FILE_TYPE_SMS
				);

				out.close();
			};

			/**
			 * @brief Writes a @ref ATEAMS::SparseMatrix to file.
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
			 * @brief Writes a @ref ATEAMS::SparseMatrix to file.
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
			 * @brief Writes a @ref ATEAMS::SparseMatrix to file.
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
