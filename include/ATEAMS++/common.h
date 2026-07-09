
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
	/** Default indexing type. */
	typedef uint32_t INDEX;

	/** Represents an element of \f$\Z/p\Z\f$. See @ref Zp. */
	typedef unsigned long FINITE;

	/** Alias of `SparseRREF::rat_t`; represents an element of \f$\Q\f$. See @ref Q. */
	typedef SparseRREF::rat_t RATIONAL;

	/** Type alias for a SparseRREF ring type. */
	typedef SparseRREF::field_t abstractRing;

	/**
	 * @class Ring
	 * @brief Generic algebraic ring.
	 * 
	 * @var Ring::ring
	 * 	Abstract (SparseRREF) ring underlying this Ring.
	 */
	struct Ring {
		public:
			abstractRing ring;
	};

	/**
	 * @class Q
	 * @brief Rationals \f$\Q\f$.
	 */
	struct Q : public Ring {
		public:
			/**
			 * @brief Data type @ref RATIONAL.
			 */
			typedef RATIONAL dtype;

			/**
			 * @brief Constructor.
			 * @param characteristic Characteristic of the field; superfluous.
			 */
			Q(int characteristic) {
				this->ring = abstractRing(SparseRREF::FIELD_QQ);
			};

			/** 
			 * @brief Constructor.
			 */
			Q() {
				this->ring = abstractRing(SparseRREF::FIELD_QQ);
			};
	};

	/**
	 * @brief Quotient ring \f$\Z/p\Z\f$.
	 */
	struct Zp : public Ring {
		public:
			/**
			 * @brief Data type @ref FINITE.
			 */
			typedef FINITE dtype;

			/**
			 * @brief Constructor.
			 * @param characteristic Characteristic (modulus) of the ring.
			 */
			Zp(int characteristic) {
				this->ring = abstractRing(SparseRREF::FIELD_Fp, characteristic);
			};
	};

	/**
	 * @brief Quotient ring \f$\Z/2\Z\f$.
	 */
	struct Z2 : public Ring {
		public:
			/**
			 * @brief Data type @ref FINITE.
			 */
			typedef FINITE dtype;

			/**
			 * @brief Constructor.
			 */
			Z2() { this->ring = abstractRing(SparseRREF::FIELD_Fp, 2); };
	};

	/**
	 * @brief Generic sparse matrix.
	 * @tparam RingLike A @ref Ring, like @ref Zp or @ref Q.
	 */
	template <typename RingLike>
	using SparseMatrix = SparseRREF::sparse_mat<typename RingLike::dtype,INDEX>;

	/**
	 * @brief Generic sparse vector.
	 * @tparam RingLike A @ref Ring, like @ref Zp or @ref Q.
	 */
	template <typename RingLike>
	using SparseVector = SparseRREF::sparse_vec<typename RingLike::dtype,INDEX>;

	/**
	 * @brief Vector of @ref SparseMatrix.
	 * @tparam RingLike A @ref Ring, like @ref Zp or @ref Q.
	 */
	template <typename RingLike>
	using SparseMatrices = std::vector<SparseMatrix<RingLike>>;

	/**
	 * @brief Vector of @ref SparseVector.
	 * @tparam RingLike A @ref Ring, like @ref Zp or @ref Q.
	 */
	template <typename RingLike>
	using SparseVectors = std::vector<SparseVector<RingLike>>;

	/**
	 * @brief Alias of @ref SparseVector.
	 * @tparam RingLike A @ref Ring, like @ref Zp or @ref Q.
	 */
	template <typename RingLike>
	using SparseBasis = std::vector<SparseVector<RingLike>>;

	/**
	 * @brief Vector of @ref SparseBasis.
	 * @tparam RingLike A @ref Ring, like @ref Zp or @ref Q.
	 */
	template <typename RingLike>
	using SparseBases = std::vector<SparseBasis<RingLike>>;

	/**
	 * @brief Alias of `std::vector`.
	 * @tparam T Numeric type.
	 */
	template <typename T>
	using DenseVector = std::vector<T,std::allocator<T>>;

	/** @cond */
	/** Flat boundary matrix. */
	typedef std::vector<std::vector<int>> FlatBoundaryMatrix;

	/** SparseRREF index of pivots. */
	typedef SparseRREF::pivot_t<INDEX> SparsePivot;

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

				SparseRREF::sparse_mat_write<Q::dtype,std::ofstream,ATEAMS::INDEX>(
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

				SparseRREF::sparse_mat_write<Zp::dtype,std::ofstream,ATEAMS::INDEX>(
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
