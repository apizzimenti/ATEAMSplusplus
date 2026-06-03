/*
	Copyright (C) 2024-2025 Zhenjie Li (Li, Zhenjie)

	This file is part of SparseRREF. The SparseRREF is free software:
	you can redistribute it and/or modify it under the terms of the MIT
	License.
*/

#ifndef SPARSE_RREF_H
#define SPARSE_RREF_H

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstring>
#include <execution>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>
#include <variant>

#include "thread_pool.hpp"

// mmap
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif

#ifdef USE_MIMALLOC
#include "mimalloc.h"
#endif

namespace SparseRREF {
	// version
	static const char version[] = "v0.3.6";
	static const int version_major = 0;
	static const int version_minor = 3;
	static const int version_patch = 6;

	enum SPARSE_FILE_TYPE {
		SPARSE_FILE_TYPE_PLAIN,
		SPARSE_FILE_TYPE_SMS,
		SPARSE_FILE_TYPE_MTX,
		SPARSE_FILE_TYPE_WXF,
		SPARSE_FILE_TYPE_BIN
	};

	// Memory management
#ifdef USE_MIMALLOC
	template <typename T>
	inline T* s_malloc(const size_t size) {
		return (T*)mi_malloc(size * sizeof(T));
	}

	template <typename T>
	inline void s_free(T* s) {
		mi_free(s);
	}

	template <typename T>
	inline T* s_realloc(T* s, const size_t size) {
		return (T*)mi_realloc(s, size * sizeof(T));
	}

	template <typename T>
	inline T* s_expand(T* s, const size_t size) {
		return (T*)mi_expand(s, size * sizeof(T));
	}
#else
	template <typename T>
	inline T* s_malloc(const size_t size) {
		return (T*)std::malloc(size * sizeof(T));
	}

	template <typename T>
	inline void s_free(T* s) {
		std::free(s);
	}

	template <typename T>
	inline T* s_realloc(T* s, const size_t size) {
		return (T*)std::realloc(s, size * sizeof(T));
	}

	template <typename T>
	inline T* s_expand(T* s, const size_t size) {
		return (T*)std::realloc(s, size * sizeof(T));
	}
#endif

	// memory info
#ifdef USE_MIMALLOC
	struct memory_info_t {
		size_t elapsed_msecs;
		size_t user_msecs;
		size_t system_msecs;
		size_t current_rss;
		size_t peak_rss;
		size_t current_commit;
		size_t peak_commit;
		size_t page_faults;
	};

	inline memory_info_t get_memory_info() {
		memory_info_t info;
		mi_process_info(
			&info.elapsed_msecs,
			&info.user_msecs,
			&info.system_msecs,
			&info.current_rss,
			&info.peak_rss,
			&info.current_commit,
			&info.peak_commit,
			&info.page_faults);
		return info;
	}

	inline void reset_memory_info() {
		mi_stats_reset();
	}
#endif

	template <typename T>
	void s_copy(T* des, const T* ini, const size_t size) {
		if (des == ini)
			return;
		std::copy(ini, ini + size, des);
	}

	// thread pool
	using thread_pool = BS::thread_pool<>;
	inline size_t thread_id() { return BS::this_thread::get_index().value(); }

	// rref_option
	// method 0: right and left search
	// method 1: only right search
	// method 2: hibrid
	struct rref_option {
		bool verbose = false;
		bool shrink_memory = false;
		std::atomic<bool> abort = false;
		bool is_back_sub = true;
		bool eliminate_one_nnz = true;
		int method = 0;
		int print_step = 100;
		std::function<int64_t(int64_t)> col_weight = [](int64_t i) { return i; };
		thread_pool pool = thread_pool(1); // default: thread pool with 1 thread
	};
	using rref_option_t = rref_option[1];

	inline size_t ctz(uint64_t x) { return std::countr_zero(x); }
	inline size_t clz(uint64_t x) { return std::countl_zero(x); }
	inline size_t popcount(uint64_t x) { return std::popcount(x); }

	template <typename T>
	inline uint8_t minimal_signed_bits(T x) noexcept {
		if (x >= INT8_MIN && x <= INT8_MAX) return 0;
		if (x >= INT16_MIN && x <= INT16_MAX) return 1;
		if (x >= INT32_MIN && x <= INT32_MAX) return 2;
		return 3; // for int64_t
	}

	template <typename T>
	inline uint8_t minimal_unsigned_bits(T x) noexcept {
		if (x <= UINT8_MAX) return 0;
		if (x <= UINT16_MAX) return 1;
		if (x <= UINT32_MAX) return 2;
		return 3; // for uint64_t
	}

	// define a special sval for index types, used for flags
	template <typename T> requires (std::is_integral_v<T>)
		constexpr T index_sval() {
		if constexpr (std::is_signed_v<T>) {
			return (T)(-1);
		}
		else {
			return std::numeric_limits<T>::max();
		}
	}

	// string
	inline void delete_space(std::string& str) {
		str.erase(std::remove_if(str.begin(), str.end(),
			[](unsigned char x) { return std::isspace(x); }),
			str.end());
	}

	inline std::vector<std::string> split_string(const std::string& s, const std::string delim) {
		size_t start = 0;
		size_t end = s.find(delim);
		std::vector<std::string> result;
		while (end != std::string::npos) {
			result.push_back(s.substr(start, end - start));
			start = end + delim.length();
			end = s.find(delim, start);
		}
		result.push_back(s.substr(start, end));
		return result;
	}

	inline uint64_t string_to_ull(std::string_view sv) {
		uint64_t result;
		auto [ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), result);
		if (ec != std::errc()) {
			throw std::runtime_error("Failed to parse number");
		}
		return result;
	}

	inline void ustr_write(const std::filesystem::path file, const std::vector<uint8_t>& str) {
		std::ofstream ofs(file, std::ios::binary);
		if (!ofs.is_open()) {
			std::cerr << "Error: ustr_write: file not open." << std::endl;
			return;
		}
		ofs.write((const char*)str.data(), str.size());
		ofs.close();
	}

	// time
	inline std::chrono::system_clock::time_point clocknow() {
		return std::chrono::system_clock::now();
	}

	inline double usedtime(std::chrono::system_clock::time_point start,
		std::chrono::system_clock::time_point end) {
		auto duration =
			std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		return ((double)duration.count() * std::chrono::microseconds::period::num /
			std::chrono::microseconds::period::den);
	}

	// some algorithms
	template <typename T> std::vector<T> difference(const std::vector<T>& l) {
		std::vector<T> result;
		for (size_t i = 1; i < l.size(); i++) {
			result.push_back(l[i] - l[i - 1]);
		}
		return result;
	}

	template <typename T>
	int lexico_compare(const T* a, const T* b, const size_t len) {
		for (size_t i = 0; i < len; i++) {
			if (a[i] < b[i])
				return -1;
			if (a[i] > b[i])
				return 1;
		}
		return 0;
	}

	template <typename T>
	int lexico_compare(const std::vector<T>& a, const std::vector<T>& b) {
		return lexico_compare(a.data(), b.data(), a.size());
	}

	template <typename T>
	int lexico_compare(const T* a, const T* b, const std::vector<size_t>& perm) {
		for (auto i : perm) {
			if (a[i] < b[i])
				return -1;
			if (a[i] > b[i])
				return 1;
		}
		return 0;
	}

	// mulit for
	template <typename Func>
	void multi_for(
		const std::vector<size_t>& start,
		const std::vector<size_t>& end,
		Func&& func) {

		if (start.size() != end.size()) {
			std::cerr << "Error: start and end size not match." << std::endl;
			exit(1);
		}

		std::vector<size_t> index(start);
		size_t nt = start.size();

		while (true) {
			func(index);

			for (int i = nt - 1; i > -2; i--) {
				if (i == -1)
					return;
				index[i]++;
				if (index[i] < end[i])
					break;
				index[i] = start[i];
			}
		}
	}

	constexpr std::array<uint64_t, 64> make_mask_table() {
		std::array<uint64_t, 64> table{};
		for (size_t i = 0; i < 64; ++i) {
			table[i] = uint64_t(1) << i;
		}
		return table;
	}

	alignas(64) constexpr auto mask_table = make_mask_table();

	// bit_array
	struct bit_array {
		std::vector<uint64_t> data;

		bit_array() {}
		~bit_array() {}

		void clear() {
			std::fill(data.begin(), data.end(), 0);
		}

		void resize(const size_t size) {
			data.resize(size / 64 + 1);
			clear();
		}

		bit_array(const size_t size) {
			resize(size);
			clear();
		}

		void insert(const size_t val) {
			auto idx = val >> 6;
			auto pos = val & 63;
			data[idx] |= mask_table[pos];
		}

		bool test(const size_t val) const {
			auto idx = val >> 6;
			auto pos = val & 63;
			return data[idx] & mask_table[pos];
		}

		void erase(const size_t val) {
			auto idx = val >> 6;
			auto pos = val & 63;
			data[idx] &= ~mask_table[pos];
		}

		void xor_insert(const size_t val) {
			auto idx = val >> 6;
			auto pos = val & 63;
			data[idx] ^= mask_table[pos];
		}

		void set(const size_t val, const bool b) {
			auto idx = val >> 6;
			auto pos = val & 63;
			data[idx] = (data[idx] & ~mask_table[pos]) | (uint64_t(b) << pos);
		}

		bool operator[](const size_t idx) const {
			return test(idx);
		}

		size_t nnz() const {
			size_t nz = 0;
			for (auto& bb : data) {
				nz += popcount(bb);
			}
			return nz;
		}

		std::vector<size_t> nonzero() const {
			auto nz = nnz();
			if (nz == 0)
				return std::vector<size_t>();

			const size_t bitset_size = 64;
			std::vector<size_t> result;
			result.reserve(nz);
			size_t tmp[64];
			size_t tmp_size = 0;
			for (size_t i = 0; i < data.size(); i++) {
				if (data[i] != 0) {
					tmp_size = 0;
					uint64_t c = data[i];

					while (c) {
						auto ctzpos = ctz(c);
						auto clzpos = bitset_size - 1 - clz(c);
						result.push_back(i * bitset_size + ctzpos);
						if (ctzpos == clzpos)
							break;
						tmp[tmp_size] = i * bitset_size + clzpos;
						tmp_size++;
						c = c ^ mask_table[ctzpos] ^ mask_table[clzpos];
					}
					for (size_t j = tmp_size; j > 0; j--) {
						result.push_back(tmp[j - 1]);
					}
				}
			}
			return result;
		}

		template <typename T>
		void nonzero(T* ptr) const {
			const size_t bitset_size = 64;
			size_t pos = 0;
			for (size_t i = 0; i < data.size(); i++) {
				if (data[i] != 0) {
					uint64_t c = data[i];
					size_t base = i * bitset_size;

					while (c) {
						auto ctzpos = ctz(c);
						ptr[pos] = base + ctzpos;
						pos++;
						c &= c - 1;
					}
				}
			}
		}

		// and clear data
		template <typename T>
		void nonzero_and_clear(T* ptr) {
			constexpr size_t bitset_size = 64;
			size_t pos = 0;
			for (size_t i = 0; i < data.size(); i++) {
				if (data[i] != 0) {
					uint64_t c = data[i];
					size_t base = i * bitset_size;

					while (c) {
						auto ctzpos = ctz(c);
						ptr[pos] = base + ctzpos;
						pos++;
						c &= c - 1;
					}
					data[i] = 0; // clear data
				}
			}
		}
	};

	template <typename T> inline T* binary_search(T* begin, T* end, T val) {
		auto ptr = std::lower_bound(begin, end, val);
		if (ptr == end || *ptr == val)
			return ptr;
		else
			return end;
	}

	template <typename T> inline T* lower_bound(T* begin, T* end, const T* val, size_t rank) {
		if (rank == 1)
			return std::lower_bound(begin, end, *val);

		size_t left = 0;
		size_t right = (end - begin) / rank;

		while (left < right) {
			size_t mid = left + (right - left) / 2;
			if (lexico_compare(begin + rank * mid, val, rank) < 0)
				left = mid + 1;
			else
				right = mid;
		}

		return begin + rank * left;
	}

	template <typename T> inline T* binary_search(T* begin, T* end, uint16_t rank, T* val) {
		auto ptr = SparseRREF::lower_bound(begin, end, rank, val);
		if (ptr == end || std::equal(ptr, ptr + rank, val))
			return ptr;
		else
			return end;
	}

	template <typename T>
	std::vector<T> perm_init(T n) {
		std::vector<T> perm(n);
		for (T i = 0; i < n; i++)
			perm[i] = i;
		return perm;
	}

	inline std::vector<size_t> perm_inverse(const std::vector<size_t>& perm) {
		size_t n = perm.size();
		std::vector<size_t> result(n);
		for (size_t i = 0; i < n; i++)
			result[perm[i]] = i;
		return result;
	}

	inline std::vector<size_t> random_perm(size_t n) {
		auto permutation = perm_init(n);

		std::random_device rd;
		std::mt19937 g(rd());

		std::shuffle(permutation.begin(), permutation.end(), g);
		return permutation;
	}

	inline bool is_identity_perm(const std::vector<size_t>& perm) {
		size_t n = perm.size();
		for (size_t i = 0; i < n; i++) {
			if (perm[i] != i)
				return false;
		}
		return true;
	}

	template <typename T>
	void permute(const std::vector<size_t>& P, T* A, size_t block_size = 1) {
		// If P is already sorted, no need to permute
		if (std::is_sorted(P.begin(), P.end())) {
			return;
		}

		std::vector<bool> visited(P.size(), false);

		auto permute_it = [&](auto& temp_block) {
			for (size_t i = 0; i < P.size(); ++i) {
				if (visited[i] || P[i] == i)
					continue;

				size_t j = i;
				for (size_t k = 0; k < block_size; ++k)
					temp_block[k] = std::move(A[j * block_size + k]);

				while (!visited[j]) {
					visited[j] = true;
					size_t k = P[j];

					if (k == i) {
						for (size_t m = 0; m < block_size; ++m)
							A[j * block_size + m] = std::move(temp_block[m]);
						break;
					}

					for (size_t m = 0; m < block_size; ++m)
						A[j * block_size + m] = std::move(A[k * block_size + m]);
					j = k;
				}
			}
			};

		if (block_size < 32) {
			// Use a fixed-size array for small block sizes
			T temp_block[32];
			permute_it(temp_block);
		}
		else {
			std::vector<T> temp_block(block_size);
			permute_it(temp_block);
		}

	}

	template <typename T>
	void permute(const std::vector<size_t>& P, std::vector<T>& A, size_t block_size = 1) {
		permute(P, A.data(), block_size);
	}

	inline std::vector<size_t> swap_perm(size_t a, size_t b, size_t n) {
		std::vector<size_t> perm(n);
		for (size_t i = 0; i < n; i++)
			perm[i] = i;
		perm[a] = b;
		perm[b] = a;
		return perm;
	}

	// some print helpers
	template <typename T, typename S>
	void print_vec(S&& os, const std::vector<T> v, const std::string_view delim = " ") {
		for (size_t i = 0; i < v.size(); i++) {
			os << v[i];
			if (i + 1 != v.size())
				os << delim;
		}
	}
	template <typename T, typename S>
	void print_vec(S&& os, const std::span<T> v, const std::string_view delim = " ") {
		for (size_t i = 0; i < v.size(); i++) {
			os << v[i];
			if (i + 1 != v.size())
				os << delim;
		}
	}

	template <typename T>
	void print_vec(T v, const std::string_view delim = " ") {
		print_vec(std::cout, v, delim);
	}

	struct MMapFile {
		std::string_view view;
#if defined(_WIN32)
		HANDLE hFile = INVALID_HANDLE_VALUE;
		HANDLE hMap = nullptr;
#else
		int fd = -1;
#endif
		size_t size = 0;

		~MMapFile() {
#if defined(_WIN32)
			if (!view.empty()) UnmapViewOfFile(view.data());
			if (hMap) CloseHandle(hMap);
			if (hFile != INVALID_HANDLE_VALUE) CloseHandle(hFile);
#else
			if (!view.empty()) munmap((void*)view.data(), size);
			if (fd >= 0) close(fd);
#endif
		}
	};

	inline bool mmap_file(const char* path, MMapFile& out) {
#if defined(_WIN32)
		out.hFile = CreateFileA(path, GENERIC_READ, FILE_SHARE_READ, nullptr,
			OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
		if (out.hFile == INVALID_HANDLE_VALUE) return false;

		LARGE_INTEGER li;
		if (!GetFileSizeEx(out.hFile, &li)) return false;
		out.size = static_cast<size_t>(li.QuadPart);

		out.hMap = CreateFileMappingA(out.hFile, nullptr, PAGE_READONLY, 0, 0, nullptr);
		if (!out.hMap) return false;

		LPVOID data = MapViewOfFile(out.hMap, FILE_MAP_READ, 0, 0, 0);
		if (!data) return false;

		out.view = std::string_view(static_cast<const char*>(data), out.size);
		return true;
#else
		out.fd = open(path, O_RDONLY);
		if (out.fd < 0) return false;

		struct stat st;
		if (fstat(out.fd, &st) < 0) return false;
		out.size = static_cast<size_t>(st.st_size);

		void* data = mmap(nullptr, out.size, PROT_READ, MAP_PRIVATE, out.fd, 0);
		if (data == MAP_FAILED) return false;

		out.view = std::string_view(static_cast<const char*>(data), out.size);
		return true;
#endif
	}

	inline std::vector<uint8_t> file_to_ustr(const std::filesystem::path filename) {
		if (!std::filesystem::exists(filename)) {
			std::cerr << "Error: File does not exist!" << std::endl;
			return std::vector<uint8_t>();
		}
		std::ifstream file(filename, std::ios::binary | std::ios::ate);
		std::streamsize size = file.tellg();
		file.seekg(0, std::ios::beg);

		std::vector<uint8_t> buffer(size);
		if (!file.read((char*)buffer.data(), size)) {
			std::cerr << "Failed to read file!" << std::endl;
			return std::vector<uint8_t>();
		}

		return buffer;
	}

} // namespace SparseRREF

#endif