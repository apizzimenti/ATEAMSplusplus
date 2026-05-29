/*
	Copyright (C) 2025 Zhenjie Li (Li, Zhenjie)

	This file is part of SparseRREF. The SparseRREF is free software:
	you can redistribute it and/or modify it under the terms of the MIT
	License.
*/


#ifndef SPARSE_TYPE_H
#define SPARSE_TYPE_H

#include "sparse_rref.h"

namespace SparseRREF {
	enum SPARSE_TYPE {
		SPARSE_CSR, // Compressed sparse row
		SPARSE_COO, // Coordinate list
		SPARSE_LR  // List of rows
	};

	template <Flint::builtin_integral index_t = int>
	struct pivot_t {
		index_t r;
		index_t c;
	};

	// sparse vector
	template <typename T, Flint::builtin_integral index_t = int> struct sparse_vec {
		index_t* indices = nullptr;
		T* entries = nullptr;
		size_t _nnz = 0;
		size_t _alloc = 0;

		struct de_iterator_ref {
			index_t& ind;
			T& val;
		};

		struct iterator {
			index_t* ind_ptr;
			T* val_ptr;

			iterator& operator++() { ind_ptr++; val_ptr++; return *this; }
			iterator operator++(int) { iterator tmp = *this; ind_ptr++; val_ptr++; return tmp; }
			iterator& operator--() { ind_ptr--; val_ptr--; return *this; }
			iterator operator--(int) { iterator tmp = *this; ind_ptr--; val_ptr--; return tmp; }
			iterator& operator+=(size_t n) { ind_ptr += n; val_ptr += n; return *this; }
			iterator& operator-=(size_t n) { ind_ptr -= n; val_ptr -= n; return *this; }
			iterator operator+(size_t n) const { iterator tmp = *this; tmp += n; return tmp; }
			iterator operator-(size_t n) const { iterator tmp = *this; tmp -= n; return tmp; }
			bool operator==(const iterator& other) const { return ind_ptr == other.ind_ptr; }
			bool operator!=(const iterator& other) const { return ind_ptr != other.ind_ptr; }

			de_iterator_ref operator*() const { return { *ind_ptr, *val_ptr }; }
		};

		// functions of iterator 
		iterator begin() { return { indices, entries }; }
		iterator end() { return { indices + _nnz, entries + _nnz }; }
		iterator begin() const { return { indices, entries }; }
		iterator end() const { return { indices + _nnz, entries + _nnz }; }
		iterator cbegin() const { return { indices, entries }; }
		iterator cend() const { return { indices + _nnz, entries + _nnz }; }

		auto index_span() const { return std::span<index_t>(indices, _nnz); }
		auto entry_span() const { return std::span<T>(entries, _nnz); }

		// C++23 is needed for zip_view
		// auto index_view() const { return std::ranges::subrange(indices, indices + _nnz); }
		// auto entry_view() const { return std::ranges::subrange(entries, entries + _nnz); }
		// auto combine_view() const { return std::ranges::zip_view(index_view(), entry_view()); }

		sparse_vec() {
			indices = nullptr;
			entries = nullptr;
			_nnz = 0;
			_alloc = 0;
		}

		void clear() {
			if (_alloc == 0)
				return;
			s_free(indices);
			indices = nullptr;
			for (size_t i = 0; i < _alloc; i++)
				entries[i].~T();
			s_free(entries);
			entries = nullptr;
			_alloc = 0;
			_nnz = 0;
		}

		~sparse_vec() {
			clear();
		}

		void reserve(size_t n, const bool is_copy = true) {
			if (n == _alloc)
				return;

			if (n == 0) {
				clear();
				return;
			}

			if (_alloc == 0) {
				indices = s_malloc<index_t>(n);
				entries = s_malloc<T>(n);
				for (size_t i = 0; i < n; i++)
					new (entries + i) T();
				_alloc = n;
				return;
			}

			// when expanding or using realloc, sometimes we need to copy the old data
			// if is_copy is false, we do not make sure that the old data is copied to 
			// the new memory, it is useful when we just want to enlarge the memory
			if (!is_copy && n > _alloc) {
				auto ii = s_expand(indices, n);
				auto ee = s_expand(entries, n);
				if (ii == nullptr) {
					s_free(indices);
					indices = s_malloc<index_t>(n);
				}
				else {
					indices = ii;
				}
				if (ee == nullptr) {
					for (size_t i = 0; i < _alloc; i++) {
						entries[i].~T();
					}
					s_free(entries);
					entries = s_malloc<T>(n);
					for (size_t i = 0; i < n; i++) {
						new (entries + i) T();
					}
				}
				else {
					entries = ee;
					for (size_t i = _alloc; i < n; i++) {
						new (entries + i) T();
					}
				}

				_alloc = n;
				_nnz = 0;
				return;
			}

			indices = s_realloc(indices, n);

			if (n < _alloc) {
				for (size_t i = n; i < _alloc; i++)
					entries[i].~T();
				entries = s_realloc<T>(entries, n);
			}
			else {
				entries = s_realloc<T>(entries, n);
				for (size_t i = _alloc; i < n; i++)
					new (entries + i) T();
			}

			_alloc = n;
		}

		inline void zero() { _nnz = 0; }
		inline void resize(size_t n) { _nnz = n; }

		inline void copy(const sparse_vec& l) {
			if (this == &l)
				return;
			if (_alloc < l._nnz)
				reserve(l._nnz);
			for (size_t i = 0; i < l._nnz; i++) {
				indices[i] = l.indices[i];
				entries[i] = l.entries[i];
			}
			_nnz = l._nnz;
		}

		inline sparse_vec(const sparse_vec& l) { copy(l); }

		inline size_t nnz() const { return _nnz; }
		inline size_t size() const { return _nnz; }
		inline size_t alloc() const { return _alloc; }

		sparse_vec(sparse_vec&& l) noexcept {
			indices = l.indices;
			entries = l.entries;
			_nnz = l._nnz;
			_alloc = l._alloc;
			l.indices = nullptr;
			l.entries = nullptr;
			l._nnz = 0;
			l._alloc = 0;
		}

		sparse_vec& operator=(const sparse_vec& l) {
			if (this == &l)
				return *this;

			copy(l);
			return *this;
		}

		sparse_vec& operator=(sparse_vec&& l) noexcept {
			if (this == &l)
				return *this;

			clear();
			indices = l.indices;
			entries = l.entries;
			_nnz = l._nnz;
			_alloc = l._alloc;
			l.indices = nullptr;
			l.entries = nullptr;
			l._nnz = 0;
			l._alloc = 0;
			return *this;
		}

		// this comparison does not clear zero elements / sort the order of indices
		bool operator==(const sparse_vec& l) const {
			if (this == &l)
				return true;
			if (_nnz != l._nnz)
				return false;
			return std::equal(indices, indices + _nnz, l.indices) 
				&& std::equal(entries, entries + _nnz, l.entries);
		}

		// this comparison will clear zero elements / sort the order of indices
		bool is_equal_to(const sparse_vec& l) const {
			if (this == &l)
				return true;
			auto this_temp = *this;
			auto other_temp = l;
			this_temp.canonicalize();
			other_temp.canonicalize();
			if (this_temp._nnz != other_temp._nnz)
				return false;
			this_temp.sort_indices();
			other_temp.sort_indices();
			return std::equal(this_temp.indices, this_temp.indices + this_temp._nnz, other_temp.indices)
				&& std::equal(this_temp.entries, this_temp.entries + this_temp._nnz, other_temp.entries);
		}

		inline void push_back(const index_t index, const T& val) {
			if (_nnz + 1 > _alloc)
				reserve((1 + _alloc) * 2); // +1 to avoid _alloc = 0
			indices[_nnz] = index;
			entries[_nnz] = val;
			_nnz++;
		}

		inline void push_back(const index_t index, T&& val) {
			if (_nnz + 1 > _alloc)
				reserve((1 + _alloc) * 2); // +1 to avoid _alloc = 0
			indices[_nnz] = index;
			entries[_nnz] = std::move(val);
			_nnz++;
		}

		index_t& operator()(const size_t pos) { return indices[pos]; }
		const index_t& operator()(const size_t pos) const { return indices[pos]; }
		T& operator[](const size_t pos) { return entries[pos]; }
		const T& operator[](const size_t pos) const { return entries[pos]; }

		// conversion functions
		template <typename U = T> requires std::is_integral_v<U> || std::is_same_v<U, int_t>
		operator sparse_vec<rat_t, index_t>() {
			sparse_vec<rat_t, index_t> result;
			result.reserve(_nnz);
			result.resize(_nnz);
			for (size_t i = 0; i < _nnz; i++) {
				result.indices[i] = indices[i];
				result.entries[i] = entries[i];
			}
			return result;
		}

		template <typename U = T> requires std::is_integral_v<U>
		operator sparse_vec<int_t, index_t>() {
			sparse_vec<int_t, index_t> result;
			result.reserve(_nnz);
			result.resize(_nnz);
			for (size_t i = 0; i < _nnz; i++) {
				result.indices[i] = indices[i];
				result.entries[i] = entries[i];
			}
			return result;
		}

		template <typename U = T> requires std::is_same_v<U, rat_t>
		sparse_vec<ulong, index_t> operator%(const nmod_t mod) const {
			sparse_vec<ulong, index_t> result;
			result.reserve(_nnz);
			result.resize(_nnz);
			for (size_t i = 0; i < _nnz; i++) {
				result.indices[i] = indices[i];
				result.entries[i] = entries[i] % mod;
			}
			return result;
		}

		template <typename U = T> requires std::is_same_v<U, rat_t>
		sparse_vec<ulong, index_t> operator%(const ulong p) const {
			sparse_vec<ulong, index_t> result;
			nmod_t mod;
			nmod_init(&mod, p);
			result.reserve(_nnz);
			result.resize(_nnz);
			for (size_t i = 0; i < _nnz; i++) {
				result.indices[i] = indices[i];
				result.entries[i] = entries[i] % mod;
			}
			return result;
		}

		void canonicalize() {
			size_t new_nnz = 0;
			for (size_t i = 0; i < _nnz; i++) {
				if (entries[i] != 0) {
					if (new_nnz != i) {
						indices[new_nnz] = indices[i];
						entries[new_nnz] = entries[i];
					}
					new_nnz++;
				}
			}
			_nnz = new_nnz;
		}

		void sort_indices() {
			if (_nnz <= 1 || std::is_sorted(indices, indices + _nnz))
				return;

			auto perm = perm_init(_nnz);
			std::sort(perm.begin(), perm.end(), [&](index_t a, index_t b) {
				return indices[a] < indices[b];
				});

			permute(perm, indices);
			permute(perm, entries);
		}

		void compress() {
			canonicalize();
			sort_indices();
			reserve(_nnz);
		}
	};

	template <Flint::builtin_integral index_t> struct sparse_vec<bool, index_t> {
		index_t* indices = nullptr;
		size_t _nnz = 0;
		size_t _alloc = 0;

		auto index_span() const { return std::span<index_t>(indices, _nnz); }

		sparse_vec() {
			indices = nullptr;
			_nnz = 0;
			_alloc = 0;
		}

		void clear() {
			if (_alloc != 0)
				s_free(indices);
			_alloc = 0;
			_nnz = 0;
		}

		~sparse_vec() { clear(); }

		void reserve(size_t n) {
			if (n == _alloc)
				return;

			if (n == 0) {
				clear();
				return;
			}

			if (_alloc == 0) {
				indices = s_malloc<index_t>(n);
				_alloc = n;
				return;
			}

			indices = s_realloc(indices, n);
			_alloc = n;
		}

		void resize(size_t n) { _nnz = n; }

		inline void copy(const sparse_vec& l) {
			if (this == &l)
				return;
			if (_alloc < l._nnz)
				reserve(l._nnz);
			for (size_t i = 0; i < l._nnz; i++) {
				indices[i] = l.indices[i];
			}
			_nnz = l._nnz;
		}

		sparse_vec(const sparse_vec& l) { copy(l); }
		size_t nnz() const { return _nnz; }

		sparse_vec(sparse_vec&& l) noexcept {
			indices = l.indices;
			_nnz = l._nnz;
			_alloc = l._alloc;
			l.indices = nullptr;
			l._nnz = 0;
			l._alloc = 0;
		}

		sparse_vec& operator=(const sparse_vec& l) {
			if (this == &l)
				return *this;

			copy(l);
			return *this;
		}

		sparse_vec& operator=(sparse_vec&& l) noexcept {
			if (this == &l)
				return *this;

			clear();
			indices = l.indices;
			_nnz = l._nnz;
			_alloc = l._alloc;
			l.indices = nullptr;
			l._nnz = 0;
			l._alloc = 0;
			return *this;
		}

		void push_back(const index_t index, const bool val = true) {
			if (_nnz + 1 > _alloc)
				reserve((1 + _alloc) * 2); // +1 to avoid _alloc = 0
			indices[_nnz] = index;
			_nnz++;
		}

		index_t& operator()(const size_t pos) { return indices[pos]; }
		const index_t& operator()(const size_t pos) const { return indices[pos]; }
		void zero() { _nnz = 0; }
		void sort_indices() { std::sort(indices, indices + _nnz); }
		void canonicalize() {}
		void compress() { sort_indices(); }
	};

	template <typename T, Flint::builtin_integral index_t = int> struct sparse_mat {
		size_t nrow = 0;
		size_t ncol = 0;
		std::vector<sparse_vec<T, index_t>> rows;

		void init(size_t r, size_t c) {
			nrow = r;
			ncol = c;
			rows = std::vector<sparse_vec<T, index_t>>(r);
		}

		sparse_mat() { nrow = 0; ncol = 0; }
		~sparse_mat() {}
		sparse_mat(size_t r, size_t c) { init(r, c); }

		sparse_vec<T, index_t>& operator[](size_t i) { return rows[i]; }
		const sparse_vec<T, index_t>& operator[](size_t i) const { return rows[i]; }

		sparse_mat(const sparse_mat& l) {
			init(l.nrow, l.ncol);
			rows = l.rows;
		}

		sparse_mat(sparse_mat&& l) noexcept {
			nrow = l.nrow;
			ncol = l.ncol;
			rows = std::move(l.rows);
			l.nrow = 0;
		}

		sparse_mat& operator=(const sparse_mat& l) {
			if (this == &l)
				return *this;
			nrow = l.nrow;
			ncol = l.ncol;
			rows = l.rows;
			return *this;
		}

		sparse_mat& operator=(sparse_mat&& l) noexcept {
			if (this == &l)
				return *this;
			nrow = l.nrow;
			ncol = l.ncol;
			rows = std::move(l.rows);
			l.nrow = 0;
			return *this;
		}

		void zero() {
			for (size_t i = 0; i < nrow; i++)
				rows[i].zero();
		}

		void clear(){
			for (size_t i = 0; i < nrow; i++) {
				rows[i].clear();
			}
			std::vector<sparse_vec<T, index_t>> tmp;
			rows.swap(tmp); // clear the vector and free memory
			nrow = 0;
			ncol = 0;
		}

		size_t nnz() const {
			size_t n = 0;
			for (size_t i = 0; i < nrow; i++)
				n += rows[i].nnz();
			return n;
		}

		size_t alloc() const {
			size_t n = 0;
			for (size_t i = 0; i < nrow; i++)
				n += rows[i].alloc();
			return n;
		}

		void compress() {
			for (size_t i = 0; i < nrow; i++) {
				rows[i].compress();
			}
		}

		void clear_zero_row() {
			size_t new_nrow = 0;
			for (size_t i = 0; i < nrow; i++) {
				if (rows[i].nnz() != 0) {
					std::swap(rows[new_nrow], rows[i]);
					new_nrow++;
				}
			}
			nrow = new_nrow;
			rows.resize(nrow);
			rows.shrink_to_fit();
		}

		sparse_mat<T, index_t> transpose() {
			sparse_mat<T, index_t> res(ncol, nrow);
			for (size_t i = 0; i < ncol; i++)
				res[i].zero();

			for (size_t i = 0; i < nrow; i++) {
				for (size_t j = 0; j < rows[i].nnz(); j++) {
					res[rows[i](j)].push_back(i, rows[i][j]);
				}
			}
			return res;
		}

		// take a span of rows
		// sparse_mat.take({start, end}) returns a sparse_mat with rows indexed in [start, end)
		sparse_mat<T, index_t> take(const std::pair<size_t, size_t>& span, thread_pool* pool = nullptr) const {
			sparse_mat<T, index_t> res(span.second - span.first, ncol);
			
			if (pool == nullptr) {
				res.rows.assign(rows.begin() + span.first, rows.begin() + span.second);
			}
			else {
				pool->detach_loop(span.first, span.second, [&](size_t i) {
					res[i - span.first] = rows[i];
				});
				pool->wait();
			}
			return res;
		}

		// append other sparse_mat to this one
		void append(const sparse_mat<T, index_t>& other, thread_pool* pool = nullptr) {
			ncol = std::max(ncol, other.ncol);

			if (pool == nullptr) {
				rows.insert(rows.end(), other.rows.begin(), other.rows.end());
				nrow += other.nrow;
				return;
			}

			rows.reserve(nrow + other.nrow);
			pool->detach_loop(0, other.nrow, [&](size_t i) {
				rows.emplace_back(other.rows[i]);
			});
			pool->wait();
			nrow += other.nrow;
		}

		void append(sparse_mat<T, index_t>&& other) {
			ncol = std::max(ncol, other.ncol);

			rows.insert(rows.end(), std::make_move_iterator(other.rows.begin()), std::make_move_iterator(other.rows.end()));
			nrow += other.nrow;
			other.clear(); // clear the other matrix
		}

		// sort rows by nnz
		void sort_rows_by_nnz() {
			std::stable_sort(rows.begin(), rows.end(), [](const sparse_vec<T, index_t>& a, const sparse_vec<T, index_t>& b) {
				return a.nnz() < b.nnz();
			});
		}

		template <typename U = T> requires std::is_same_v<U, rat_t>
		sparse_mat<ulong> operator%(const nmod_t mod) const {
			sparse_mat<ulong> result(nrow, ncol);
			for (size_t i = 0; i < nrow; i++) {
				result[i] = rows[i] % mod;
			}
			return result;
		}

		template <typename U = T> requires std::is_same_v<U, rat_t>
		ulong height_bits() const {
			ulong max_height = 0;
			for (size_t i = 0; i < nrow; i++) {
				for (size_t j = 0; j < rows[i].nnz(); j++) {
					auto rr = rows[i][j].height_bits();
					if (rr > max_height)
						max_height = rr;
				}
			}

			return max_height;
		}

		// denominator bits
		template <typename U = T> requires std::is_same_v<U, rat_t>
		ulong den_bits() const {
			int_t den = 1;
			for (size_t i = 0; i < nrow; i++) {
				for (size_t j = 0; j < rows[i].nnz(); j++) {
					if (rows[i][j] != 0 && !(rows[i][j].is_den_one()))
						den = LCM(den, rows[i][j].den());
				}
			}

			return den.bits();
		}
	};

	// CSR format for sparse tensor
	template <typename T, typename index_t> struct sparse_tensor_struct {
		size_t rank;
		size_t alloc;
		index_t* colptr;
		T* valptr;
		std::vector<size_t> dims;
		std::vector<size_t> rowptr;

		using index_v = std::vector<index_t>;
		using index_p = index_t*;
		using const_index_p = const index_t*;

		//empty constructor
		sparse_tensor_struct() {
			rank = 0;
			alloc = 0;
			colptr = nullptr;
			valptr = nullptr;
		}

		// Constructor with dimensions
		// we require that rank >= 2
		void init(const std::vector<size_t>& l, size_t aoc = 8) {
			dims = l;
			rank = l.size();
			rowptr = std::vector<size_t>(l[0] + 1, 0);
			alloc = aoc;
			colptr = s_malloc<index_t>((rank - 1) * alloc);
			valptr = s_malloc<T>(alloc);
			for (size_t i = 0; i < alloc; i++)
				new (valptr + i) T();
		}

		sparse_tensor_struct(const std::vector<size_t>& l, size_t aoc = 8) {
			init(l, aoc);
		}

		// Copy constructor
		sparse_tensor_struct(const sparse_tensor_struct& l) {
			init(l.dims, l.alloc);
			std::copy(l.rowptr.begin(), l.rowptr.end(), rowptr.begin());
			std::copy(l.colptr, l.colptr + alloc * (rank - 1), colptr);
			for (size_t i = 0; i < alloc; i++)
				valptr[i] = l.valptr[i];
		}

		// Move constructor
		sparse_tensor_struct(sparse_tensor_struct&& l) noexcept {
			dims = l.dims;
			rank = l.rank;
			rowptr = l.rowptr;
			alloc = l.alloc;
			colptr = l.colptr;
			l.colptr = nullptr;
			valptr = l.valptr;
			l.valptr = nullptr;
			l.alloc = 0; // important for no repeating clear
		}

		void clear() {
			if (alloc == 0)
				return;
			for (size_t i = 0; i < alloc; i++)
				valptr[i].~T();
			s_free(valptr);
			s_free(colptr);
			valptr = nullptr;
			colptr = nullptr;
			alloc = 0;
		}

		~sparse_tensor_struct() {
			clear();
		}

		void reserve(size_t size) {
			if (size == alloc)
				return;
			if (size == 0) {
				clear();
				return;
			}
			if (alloc == 0) {
				alloc = size;
				colptr = s_malloc<index_t>(size * (rank - 1));
				valptr = s_malloc<T>(size);
				for (size_t i = 0; i < size; i++)
					new (valptr + i) T();
				return;
			}
			colptr = s_realloc<index_t>(colptr, size * (rank - 1));
			if (size > alloc) {
				valptr = s_realloc<T>(valptr, size);
				for (size_t i = alloc; i < size; i++)
					new (valptr + i) T();
			}
			else if (size < alloc) {
				for (size_t i = size; i < alloc; i++)
					valptr[i].~T();
				valptr = s_realloc<T>(valptr, size);
			}
			alloc = size;
		}

		void zero() {
			if (rank != 0)
				std::fill(rowptr.begin(), rowptr.end(), 0);
		}

		inline size_t nnz() {
			return rowptr[dims[0]];
		}

		// Copy assignment
		sparse_tensor_struct& operator=(const sparse_tensor_struct& l) {
			if (this == &l)
				return *this;
			auto nz = l.nnz();
			if (alloc == 0) {
				init(l.dims, nz);
				std::copy(l.rowptr.begin(), l.rowptr.end(), rowptr.begin());
				std::copy(l.colptr, l.colptr + nz * (rank - 1), colptr);
				std::copy(l.valptr, l.valptr + nz, valptr);
				return *this;
			}
			dims = l.dims;
			rank = l.rank;
			rowptr = l.rowptr;
			if (alloc < nz)
				reserve(nz);
			std::copy(l.colptr, l.colptr + nz * (rank - 1), colptr);
			std::copy(l.valptr, l.valptr + nz, valptr);
			return *this;
		}

		// Move assignment
		sparse_tensor_struct& operator=(sparse_tensor_struct&& l) noexcept {
			if (this == &l)
				return *this;
			clear();
			dims = l.dims;
			rank = l.rank;
			rowptr = l.rowptr;
			alloc = l.alloc;
			colptr = l.colptr;
			l.colptr = nullptr;
			valptr = l.valptr;
			l.valptr = nullptr;
			l.alloc = 0; // important for no repeating clear
			return *this;
		}

		std::vector<size_t> row_nums() {
			return SparseRREF::difference(rowptr);
		}

		size_t row_nnz(size_t i) {
			return rowptr[i + 1] - rowptr[i];
		}

		// remove zero entries, double pointer
		void canonicalize() {
			size_t nnz_now = nnz();
			size_t index = 0;
			std::vector<size_t> newrowptr(dims[0] + 1);
			newrowptr[0] = 0;
			for (size_t i = 0; i < dims[0]; i++) {
				for (size_t j = rowptr[i]; j < rowptr[i + 1]; j++) {
					if (valptr[j] != 0) {
						s_copy(colptr + index * (rank - 1), colptr + j * (rank - 1), rank - 1);
						valptr[index] = valptr[j];
						index++;
					}
				}
				newrowptr[i + 1] = index;
			}
			rowptr = newrowptr;
		}

		std::pair<index_p, T*> row(size_t i) {
			return std::make_pair(colptr + rowptr[i] * (rank - 1), valptr + rowptr[i]);
		}

		index_p entry_lower_bound(const_index_p l) {
			auto begin = row(l[0]).first;
			auto end = row(l[0] + 1).first;
			if (begin == end)
				return end;
			return SparseRREF::lower_bound(begin, end, l + 1, rank - 1);
		}

		index_p entry_lower_bound(const index_v& l) {
			return entry_lower_bound(l.data());
		}

		index_p entry_ptr(index_p l) {
			auto ptr = entry_lower_bound(l);
			auto end = row(l[0] + 1).first;
			if (ptr == end || std::equal(ptr, ptr + rank - 1, l + 1))
				return ptr;
			else
				return end;
		}

		index_p entry_ptr(const index_v& l) {
			return entry_ptr(l.data());
		}

		// unordered, push back on the end of the row
		void push_back(const index_v& l, const T& val) {
			index_t row = l[0];
			size_t nnz = this->nnz();
			if (nnz + 1 > alloc)
				reserve((alloc + 1) * 2);
			size_t index = rowptr[row + 1];
			for (size_t i = nnz; i > index; i--) {
				auto tmpptr = colptr + (i - 1) * (rank - 1);
				std::copy_backward(tmpptr, tmpptr + (rank - 1), tmpptr + 2 * (rank - 1));
				valptr[i] = valptr[i - 1];
			}
			for (size_t i = 0; i < rank - 1; i++)
				colptr[index * (rank - 1) + i] = l[i + 1];
			valptr[index] = val;
			for (size_t i = row + 1; i <= dims[0]; i++)
				rowptr[i]++;
		}

		// ordered insert
		// mode = false: insert anyway
		// mode = true: insert and replace if exist
		void insert(const index_v& l, const T& val, bool mode = true) {
			size_t trow = l[0];
			size_t nnz = this->nnz();
			if (nnz + 1 > alloc)
				reserve((alloc + 1) * 2);
			auto ptr = entry_lower_bound(l);
			size_t index = (ptr - colptr) / (rank - 1);
			bool exist = (ptr != row(trow + 1).first && std::equal(ptr, ptr + rank - 1, l.data() + 1));
			if (!exist || !mode) {
				for (size_t i = nnz; i > index; i--) {
					auto tmpptr = colptr + (i - 1) * (rank - 1);
					std::copy_backward(tmpptr, tmpptr + (rank - 1), tmpptr + 2 * (rank - 1));
					valptr[i] = valptr[i - 1];
				}
				std::copy(l.begin() + 1, l.begin() + rank, colptr + index * (rank - 1));
				valptr[index] = val;
				for (size_t i = trow + 1; i <= dims[0]; i++)
					rowptr[i]++;
				return;
			}
			valptr[index] = val;
		}

		// ordered add one value
		void insert_add(const index_v& l, const T& val) {
			size_t trow = l[0];
			size_t nnz = this->nnz();
			if (nnz + 1 > alloc)
				reserve((alloc + 1) * 2);
			auto ptr = entry_lower_bound(l);
			size_t index = (ptr - colptr) / (rank - 1);
			bool exist = (ptr != row(trow + 1).first && std::equal(ptr, ptr + rank - 1, l.data() + 1));
			if (!exist) {
				for (size_t i = nnz; i > index; i--) {
					auto tmpptr = colptr + (i - 1) * (rank - 1);
					std::copy_backward(tmpptr, tmpptr + (rank - 1), tmpptr + 2 * (rank - 1));
					valptr[i] = valptr[i - 1];
				}
				std::copy(l.begin() + 1, l.begin() + rank, colptr + index * (rank - 1));
				valptr[index] = val;
				for (size_t i = trow + 1; i <= dims[0]; i++)
					rowptr[i]++;
				return;
			}
			valptr[index] += val;
		}

		sparse_tensor_struct<T, index_t> transpose(const std::vector<size_t>& perm) {
			std::vector<size_t> l(rank);
			std::vector<size_t> lperm(rank);
			std::vector<index_t> old_coord(rank);
			std::vector<index_t> new_coord(rank);
			for (size_t i = 0; i < rank; i++)
				lperm[i] = dims[perm[i]];
			sparse_tensor_struct<T, index_t> B(lperm, nnz());
			for (size_t i = 0; i < dims[0]; i++) {
				for (size_t j = rowptr[i]; j < rowptr[i + 1]; j++) {
					old_coord[0] = i;
					auto tmpptr = colptr + j * (rank - 1);
					for (size_t k = 1; k < rank; k++)
						old_coord[k] = tmpptr[k - 1];
					for (size_t k = 0; k < rank; k++)
						new_coord[k] = old_coord[perm[k]];
					B.push_back(new_coord, valptr[j]);
				}
			}
			return B;
		}

		// multithread version use more memory and it will compress the tensor
		void sort_indices(thread_pool* pool = nullptr) {
			if (pool == nullptr) {
				for (size_t i = 0; i < dims[0]; i++) {
					size_t rownnz = rowptr[i + 1] - rowptr[i];
					std::vector<size_t> perm(rownnz);
					for (size_t j = 0; j < rownnz; j++)
						perm[j] = j;
					std::sort(std::execution::par, perm.begin(), perm.end(), [&](size_t a, size_t b) {
						auto ptra = colptr + (rowptr[i] + a) * (rank - 1);
						auto ptrb = colptr + (rowptr[i] + b) * (rank - 1);
						return lexico_compare(ptra, ptrb, rank - 1) < 0;
						});

					permute(perm, colptr + rowptr[i] * (rank - 1), rank - 1);
					permute(perm, valptr + rowptr[i]);
				}
				return;
			}

			auto nz = nnz();
			auto n_colptr = s_malloc<index_t>(nz * (rank - 1));
			auto n_valptr = s_malloc<T>(nz);

			for (size_t i = 0; i < dims[0]; i++) {
				size_t rownnz = rowptr[i + 1] - rowptr[i];
				std::vector<size_t> perm(rownnz);
				for (size_t j = 0; j < rownnz; j++)
					perm[j] = j;
				std::sort(std::execution::par, perm.begin(), perm.end(), [&](size_t a, size_t b) {
					auto ptra = colptr + (rowptr[i] + a) * (rank - 1);
					auto ptrb = colptr + (rowptr[i] + b) * (rank - 1);
					return lexico_compare(ptra, ptrb, rank - 1) < 0;
					});

				pool->detach_loop(0, rownnz, [&](size_t j) {
					auto oldptr = colptr + (rowptr[i] + perm[j]) * (rank - 1);
					auto newptr = n_colptr + (rowptr[i] + j) * (rank - 1);
					std::copy(oldptr, oldptr + (rank - 1), newptr);
					n_valptr[rowptr[i] + j] = valptr[rowptr[i] + perm[j]];
					});
				pool->wait();
			}
			s_free(colptr);
			s_free(valptr);
			colptr = n_colptr;
			valptr = n_valptr;
			alloc = nz;
		}
	};

	// define the default sparse tensor
	template <typename T, typename index_t = int, SPARSE_TYPE Type = SPARSE_COO> struct sparse_tensor;

	template <typename T, typename index_t> struct sparse_tensor<T, index_t, SPARSE_CSR> {
		sparse_tensor_struct<T, index_t> data;

		using index_v = std::vector<index_t>;
		using index_p = index_t*;

		void clear() { data.clear(); }

		sparse_tensor() {}
		~sparse_tensor() {}
		sparse_tensor(std::vector<size_t> l, size_t aoc = 8) : data(l, aoc) {}
		sparse_tensor(const sparse_tensor& l) : data(l.data) {}
		sparse_tensor(sparse_tensor&& l) noexcept : data(std::move(l.data)) {}
		sparse_tensor& operator=(const sparse_tensor& l) { data = l.data; return *this; }
		sparse_tensor& operator=(sparse_tensor&& l) noexcept { data = std::move(l.data); return *this; }

		inline size_t alloc() const { return data.alloc; }
		inline size_t rank() const { return data.rank; }
		inline size_t nnz() const { return data.rowptr[data.dims[0]]; }
		inline auto& rowptr() const { return data.rowptr; }
		inline auto& colptr() const { return data.colptr; }
		inline auto& valptr() const { return data.valptr; }
		inline auto& dims() const { return data.dims; }
		inline size_t dim(size_t i) const { return data.dims[i]; }
		index_p index(size_t i) const { return data.colptr + i * (rank() - 1); }
		T& val(size_t i) const { return data.valptr[i]; }
		inline void zero() { data.zero(); }
		inline void insert(const index_v& l, const T& val, bool mode = true) { data.insert(l, val, mode); }
		inline void push_back(const index_v& l, const T& val) { data.push_back(l, val); }
		inline void canonicalize() { data.canonicalize(); }
		inline void sort_indices(thread_pool* pool = nullptr) { data.sort_indices(pool); }
		inline void reserve(size_t size) { data.reserve(size); }
		inline sparse_tensor transpose(const std::vector<size_t>& perm) {
			sparse_tensor B;
			B.data = data.transpose(perm);
			return B;
		}

		void convert_from_COO(const sparse_tensor<T, index_t, SPARSE_COO>& l) {
			data.dims = l.data.dims;
			data.rank = l.data.rank;
			auto nnz = l.nnz();
			if (alloc() < nnz)
				reserve(nnz);
			std::copy(l.data.valptr, l.data.valptr + nnz, data.valptr);
			auto newrank = data.rank - 1;
			for (size_t i = 0; i < newrank; i++)
				data.dims[i] = data.dims[i + 1];
			data.dims.resize(newrank);
			// then recompute the rowptr and colptr
			// first compute nnz for each row
			data.rowptr.resize(data.dims[0] + 1, 0);
			for (size_t i = 0; i < nnz; i++) {
				auto oldptr = l.data.colptr + i * newrank;
				auto nowptr = data.colptr + i * (newrank - 1);
				data.rowptr[oldptr[0] + 1]++;
				for (size_t j = 0; j < newrank - 1; j++)
					nowptr[j] = oldptr[j + 1];
			}
			for (size_t i = 0; i < data.dims[0]; i++)
				data.rowptr[i + 1] += data.rowptr[i];
			data.rank = newrank;
		}

		void move_from_COO(sparse_tensor<T, index_t, SPARSE_COO>&& l) noexcept {
			// use move to avoid memory problems, then no need to mention l
			data = std::move(l.data);
			// remove the first dimension
			auto newrank = data.rank - 1;
			for (size_t i = 0; i < newrank; i++)
				data.dims[i] = data.dims[i + 1];
			data.dims.resize(newrank);
			// then recompute the rowptr and colptr
			// first compute nnz for each row
			std::vector<size_t> rowptr(data.dims[0] + 1, 0);
			auto nnz = data.rowptr[1];
			for (size_t i = 0; i < nnz; i++) {
				auto oldptr = data.colptr + i * newrank;
				auto nowptr = data.colptr + i * (newrank - 1);
				rowptr[oldptr[0] + 1]++;
				for (size_t j = 0; j < newrank - 1; j++)
					nowptr[j] = oldptr[j + 1];
			}
			for (size_t i = 0; i < data.dims[0]; i++)
				rowptr[i + 1] += rowptr[i];
			data.rowptr = rowptr;
			data.colptr = s_realloc<index_t>(data.colptr, nnz * (newrank - 1));
			data.rank = newrank;
		}

		// constructor from COO
		sparse_tensor(const sparse_tensor<T, index_t, SPARSE_COO>& l) { convert_from_COO(l); }
		sparse_tensor& operator=(const sparse_tensor<T, index_t, SPARSE_COO>& l) {
			convert_from_COO(l);
			return *this;
		}

		// suppose that COO tensor is sorted
		sparse_tensor(sparse_tensor<T, index_t, SPARSE_COO>&& l) noexcept {
			move_from_COO(std::move(l));
		}

		sparse_tensor& operator=(sparse_tensor<T, index_t, SPARSE_COO>&& l) noexcept {
			move_from_COO(std::move(l));
			return *this;
		}

		// convert a sparse matrix to a sparse tensor
		void convert_from_sparse_mat(const sparse_mat<T, index_t>& mat, thread_pool* pool = nullptr) {
			data.dims = { mat.nrow, mat.ncol };
			data.rank = 2;

			auto nnz = mat.nnz();
			if (alloc() < nnz) 
				reserve(nnz);

			// compute the rowptr
			data.rowptr.resize(mat.nrow + 1, 0);
			for (size_t i = 0; i < mat.nrow; i++) {
				data.rowptr[i + 1] = data.rowptr[i] + mat[i].nnz();
			}
			// copy the values and column indices
			if (pool == nullptr) {
				for (size_t i = 0; i < mat.nrow; i++) {
					std::copy(mat[i].indices, mat[i].indices + mat[i].nnz(), data.colptr + data.rowptr[i]);
					std::copy(mat[i].entries, mat[i].entries + mat[i].nnz(), data.valptr + data.rowptr[i]);
				}
			} 
			else {
				pool->detach_loop(0, mat.nrow, [&](size_t i) {
					std::copy(mat[i].indices, mat[i].indices + mat[i].nnz(), data.colptr + data.rowptr[i]);
					std::copy(mat[i].entries, mat[i].entries + mat[i].nnz(), data.valptr + data.rowptr[i]);
				});
				pool->wait();
			}
		}

		sparse_mat<T, index_t> to_sparse_mat(thread_pool* pool = nullptr) const {
			if (rank() != 2) {
				std::cerr << "sparse_tensor.to_sparse_mat: rank must be 2" << std::endl;
				return sparse_mat<T, index_t>();
			}

			sparse_mat<T, index_t> mat(data.dims[0], data.dims[1]);
			if (pool == nullptr) {
				for (size_t i = 0; i < data.dims[0]; i++) {
					auto nz = data.rowptr[i + 1] - data.rowptr[i];
					mat[i].reserve(nz);
					mat[i].resize(nz);
					std::copy(data.colptr + data.rowptr[i], data.colptr + data.rowptr[i + 1], mat[i].indices);
					std::copy(data.valptr + data.rowptr[i], data.valptr + data.rowptr[i + 1], mat[i].entries);
				}
			} 
			else {
				pool->detach_loop(0, data.dims[0], [&](size_t i) {
					auto nz = data.rowptr[i + 1] - data.rowptr[i];
					mat[i].reserve(nz);
					mat[i].resize(nz);
					std::copy(data.colptr + data.rowptr[i], data.colptr + data.rowptr[i + 1], mat[i].indices);
					std::copy(data.valptr + data.rowptr[i], data.valptr + data.rowptr[i + 1], mat[i].entries);
				});
				pool->wait();
			}
			return mat;
		}

		sparse_tensor(const sparse_mat<T, index_t>& mat, thread_pool* pool = nullptr) {
			convert_from_sparse_mat(mat, pool);
		}

		sparse_tensor& operator=(const sparse_mat<T, index_t>& mat) {
			convert_from_sparse_mat(mat);
			return *this;
		}

		// only for test
		void print_test() {
			for (size_t i = 0; i < data.dims[0]; i++) {
				for (size_t j = data.rowptr[i]; j < data.rowptr[i + 1]; j++) {
					std::cout << i << " ";
					for (size_t k = 0; k < data.rank - 1; k++)
						std::cout << (size_t)data.colptr[j * (data.rank - 1) + k] << " ";
					std::cout << " : " << data.valptr[j] << std::endl;
				}
			}
		}
	};

	template <typename index_t, typename T> struct sparse_tensor<T, index_t, SPARSE_COO> {
		sparse_tensor_struct<T, index_t> data;

		using index_v = std::vector<index_t>;
		using index_p = index_t*;
		using const_index_p = const index_t*;

		template <typename S, typename U = S> requires std::convertible_to<U, S>
		std::vector<S> prepend_num(const std::vector<S>& l, U num = 0) {
			std::vector<S> lp;
			lp.reserve(l.size() + 1);
			lp.push_back(static_cast<S>(num));
			lp.insert(lp.end(), l.begin(), l.end());
			return lp;
		}

		void clear() { data.clear(); }
		void init(const std::vector<size_t>& l, size_t aoc = 8) {
			data.init(prepend_num(l, (size_t)1), aoc);
		}

		sparse_tensor() {}
		~sparse_tensor() {}
		sparse_tensor(const std::vector<size_t>& l, size_t aoc = 8) : data(prepend_num(l, (size_t)1), aoc) {}
		sparse_tensor(const sparse_tensor& l) : data(l.data) {}
		sparse_tensor(sparse_tensor&& l) noexcept : data(std::move(l.data)) {}
		sparse_tensor& operator=(const sparse_tensor& l) { data = l.data; return *this; }
		sparse_tensor& operator=(sparse_tensor&& l) noexcept { data = std::move(l.data); return *this; }

		// for the i-th column, return the indices
		index_p index(size_t i) const { return data.colptr + i * rank(); }
		T& val(size_t i) const { return data.valptr[i]; }

		index_v index_vector(size_t i) const {
			index_v result(rank());
			for (size_t j = 0; j < rank(); j++)
				result[j] = index(i)[j];
			return result;
		}

		inline size_t alloc() const { return data.alloc; }
		inline size_t nnz() const { return data.rowptr[1]; }
		inline size_t rank() const { return data.rank - 1; }
		inline std::vector<size_t> dims() const {
			std::vector<size_t> result(data.dims.begin() + 1, data.dims.end());
			return result;
		}
		inline size_t dim(size_t i) const { return data.dims[i + 1]; }
		inline void zero() { data.zero(); }
		inline void reserve(size_t size) { data.reserve(size); }
		inline void resize(size_t new_nnz) {
			if (new_nnz > alloc())
				reserve(new_nnz);
			data.rowptr[1] = new_nnz;
		}

		// we assume that the tensor is sorted
		inline std::vector<size_t> rowptr() const {
			std::vector<size_t> result(dim(0) + 1);
			result[0] = 0;
			for (auto i = 0; i < nnz(); i++) {
				result[index(i)[0] + 1]++;
			}
			for (size_t i = 0; i < dim(0); i++)
				result[i + 1] += result[i];
			return result;
		}

		// change the dimensions of the tensor
		// it is dangerous, only for internal use
		inline void change_dims(const std::vector<size_t>& new_dims) {
			auto dims = prepend_num(new_dims, (size_t)1);
			data.dims = dims;
			data.rank = dims.size();
			data.colptr = s_realloc<index_t>(data.colptr, new_dims.size() * alloc());
		}

		inline void flatten(const std::vector<std::vector<size_t>>& pos) {
			auto r = rank();
			auto nr = pos.size();
			std::vector<index_t> newindex(nr);
			std::vector<size_t> new_dims(nr);
			auto old_dim = dims();
			auto init_ptr = data.colptr;

			// first compute new dimensions
			for (size_t i = 0; i < nr; i++) {
				new_dims[i] = 1;
				for (auto j : pos[i])
					new_dims[i] *= old_dim[j];
			}
			new_dims = prepend_num(new_dims, (size_t)1);

			for (size_t i = 0; i < nnz(); i++) {
				auto ptr = index(i);
				for (size_t j = 0; j < nr; j++) {
					newindex[j] = 0;
					for (auto k : pos[j])
						newindex[j] = newindex[j] * old_dim[k] + ptr[k];
				}
				for (size_t j = 0; j < nr; j++)
					init_ptr[i * nr + j] = newindex[j];
			}
			data.colptr = s_realloc(data.colptr, nr * nnz());

			// change the dimensions
			data.dims = new_dims;
			data.rank = nr + 1;
		}

		// reshape, for example {2,100} to {2,5,20}
		// TODO: check more examples
		inline void reshape(const std::vector<size_t>& new_dims) {
			auto old_dims = dims();
			index_t* newcolptr = s_malloc<index_t>(nnz() * new_dims.size());
			auto r = rank();

			int_t flatten_index = 0;
			int_t tmp;
			for (size_t i = 0; i < nnz(); i++) {
				auto ptr = index(i);
				flatten_index = 0;
				for (size_t j = 0; j < r; j++) {
					flatten_index *= old_dims[j];
					flatten_index += ptr[j];
				}
				for (auto j = new_dims.size(); j > 0; j--) {
					tmp = flatten_index % new_dims[j - 1];
					flatten_index /= new_dims[j - 1];
					newcolptr[i * new_dims.size() + j - 1] = tmp.to_si();
				}
			}
			s_free(data.colptr);
			data.colptr = newcolptr;
			data.dims = prepend_num(new_dims, (size_t)1);
			data.rank = new_dims.size() + 1;
		}

		inline void insert(const index_v& l, const T& val, bool mode = true) { data.insert(prepend_num(l), val, mode); }
		inline void insert_add(const index_v& l, const T& val) { data.insert_add(prepend_num(l), val); }
		void push_back(const_index_p l, const T& new_val) {
			auto n_nnz = nnz();
			if (n_nnz + 1 > data.alloc)
				reserve((data.alloc + 1) * 2);
			s_copy(index(n_nnz), l, rank());
			val(n_nnz) = new_val;
			data.rowptr[1]++; // increase the nnz
		}
		void push_back(const index_v& l, const T& new_val) { push_back(l.data(), new_val); }
		inline void canonicalize() { data.canonicalize(); }
		inline void sort_indices(thread_pool* pool = nullptr) { data.sort_indices(pool); }
		inline sparse_tensor transpose(const std::vector<size_t>& perm) {
			std::vector<size_t> perm_new(perm);
			for (auto& a : perm_new) { a++; }
			perm_new = prepend_num(perm_new, (size_t)0);
			sparse_tensor B;
			B.data = data.transpose(perm_new);
			B.sort_indices();
			return B;
		}

		sparse_mat<T, index_t> to_sparse_mat(const bool sort_ind = true, thread_pool* pool = nullptr) {
			if (rank() != 2) {
				std::cerr << "sparse_tensor.to_sparse_mat: rank must be 2" << std::endl;
				return sparse_mat<T, index_t>();
			}
			if (sort_ind)
				sort_indices(pool);
			auto r = dim(0);
			auto c = dim(1);
			auto rptr = rowptr();

			sparse_mat<T, index_t> mat(r, c);
			if (pool == nullptr) {
				for (size_t i = 0; i < r; i++) {
					auto nz = rptr[i + 1] - rptr[i];
					mat[i].reserve(nz);
					mat[i].resize(nz);
					std::copy(data.valptr + rptr[i], data.valptr + rptr[i + 1], mat[i].entries);
					// skip the first index, which is the row index
					auto ptr = index(rptr[i]) + 1;
					for (size_t j = 0; j < nz; j++)
						mat[i].indices[j] = ptr[2 * j];
				}
			}
			else {
				pool->detach_loop(0, r, [&](size_t i) {
					auto nz = rptr[i + 1] - rptr[i];
					mat[i].reserve(nz);
					mat[i].resize(nz);
					std::copy(data.valptr + rptr[i], data.valptr + rptr[i + 1], mat[i].entries);
					// skip the first index, which is the row index
					auto ptr = index(rptr[i]) + 1;
					for (size_t j = 0; j < nz; j++)
						mat[i].indices[j] = ptr[2 * j];
				});
				pool->wait();
			}
			return mat;
		}

		std::vector<size_t> gen_perm() const {
			std::vector<size_t> perm = perm_init(nnz());

			auto r = rank();
			std::sort(std::execution::par, perm.begin(), perm.end(), [&](size_t a, size_t b) {
				return lexico_compare(index(a), index(b), r) < 0;
				});
			return perm;
		}

		std::vector<size_t> gen_perm(const std::vector<size_t>& index_perm) const {
			if (index_perm.size() != rank()) {
				std::cerr << "Error: gen_perm: index_perm size is not equal to rank" << std::endl;
				exit(1);
			}

			if (std::is_sorted(index_perm.begin(), index_perm.end())) 
				return gen_perm();

			std::vector<size_t> perm = perm_init(nnz());
			std::sort(std::execution::par, perm.begin(), perm.end(), [&](size_t a, size_t b) {
				return lexico_compare(index(a), index(b), index_perm) < 0;
				});

			return perm;
		}

		void transpose_replace(const std::vector<size_t>& perm, thread_pool* pool = nullptr) {
			std::vector<size_t> new_dims(rank() + 1);
			new_dims[0] = data.dims[0];

			for (size_t i = 0; i < rank(); i++)
				new_dims[i + 1] = data.dims[perm[i] + 1];
			data.dims = new_dims;

			auto method = [&](size_t ss, size_t ee) {
				std::vector<size_t> index_new(rank());
				for (size_t i = ss; i < ee; i++) {
					auto ptr = index(i);
					for (size_t j = 0; j < rank(); j++)
						index_new[j] = ptr[perm[j]];
					std::copy(index_new.begin(), index_new.end(), ptr);
				}
				};

			if (pool == nullptr) 
				method(0, nnz());
			else {
				pool->detach_blocks(0, nnz(), method);
				pool->wait();
			}
		}

		sparse_tensor<T, index_t, SPARSE_COO> chop(size_t pos, size_t aa) const {
			std::vector<size_t> dims_new = dims();
			dims_new.erase(dims_new.begin() + pos);
			sparse_tensor<T, index_t, SPARSE_COO> result(dims_new);
			index_v index_new;
			index_new.reserve(rank() - 1);
			for (size_t i = 0; i < nnz(); i++) {
				if (index(i)[pos] != aa)
					continue;
				for (size_t j = 0; j < rank(); j++) {
					if (j != pos)
						index_new.push_back(index(i)[j]);
				}
				result.push_back(index_new, val(i));
				index_new.clear();
			}
			return result;
		}

		// constructor from CSR
		sparse_tensor(const sparse_tensor<T, index_t, SPARSE_CSR>& l) {
			data.init(prepend_num(l.dims(), (size_t)1), l.nnz());
			resize(l.nnz());

			auto r = rank();
			auto n_row = dim(0);

			// first copy the data
			s_copy(data.valptr, l.data.valptr, l.nnz());

			// then deal with the indices
			for (size_t i = 0; i < n_row; i++) {
				for (size_t j = l.data.rowptr[i]; j < l.data.rowptr[i + 1]; j++) {
					auto tmp_index = index(j);
					tmp_index[0] = i;
					for (size_t k = 0; k < r - 1; k++)
						tmp_index[k + 1] = l.data.colptr[j * (r - 1) + k];
				}
			}
		}

		sparse_tensor& operator=(const sparse_tensor<T, index_t, SPARSE_CSR>& l) {
			if (alloc() == 0) {
				init(l.dims(), l.nnz());
			}
			else {
				change_dims(l.dims());
				reserve(l.nnz());
			}

			auto r = rank();
			auto n_row = dim(0);

			// first copy the data
			s_copy(data.valptr, l.data.valptr, l.nnz());

			// then deal with the indices
			for (size_t i = 0; i < n_row; i++) {
				for (size_t j = l.data.rowptr[i]; j < l.data.rowptr[i + 1]; j++) {
					auto tmp_index = index(j);
					tmp_index[0] = i;
					for (size_t k = 0; k < r - 1; k++)
						tmp_index[k + 1] = l.data.colptr[j * (r - 1) + k];
				}
			}

			return *this;
		}

		sparse_tensor& operator=(sparse_tensor<T, index_t, SPARSE_CSR>&& l) noexcept {
			data = std::move(l.data);
			if (data.alloc == 0)
				return *this;

			auto r = data.rank;
			auto n_row = data.dims[0];

			// recompute the index
			index_t* newcolptr = s_malloc<index_t>(data.alloc * r);
			auto newcolptr_j = newcolptr;
			auto nowcolptr_j = data.colptr;
			for (size_t i = 0; i < n_row; i++) {
				for (size_t j = data.rowptr[i]; j < data.rowptr[i + 1]; j++) {
					newcolptr_j[0] = i;
					s_copy(newcolptr_j + 1, nowcolptr_j, r - 1);
					newcolptr_j += r;
					nowcolptr_j += r - 1;
				}
			}
			s_free(data.colptr);
			data.colptr = newcolptr;

			data.rowptr = { 0, data.rowptr.back() };
			data.dims = prepend_num(data.dims, (size_t)1);
			data.rank++;

			return *this;
		}

		sparse_tensor(sparse_tensor<T, index_t, SPARSE_CSR>&& l) noexcept {
			data = std::move(l.data);
			if (data.alloc == 0)
				return;

			auto r = data.rank;
			auto n_row = data.dims[0];

			// recompute the index
			index_t* newcolptr = s_malloc<index_t>(data.alloc * r);
			auto newcolptr_j = newcolptr;
			auto nowcolptr_j = data.colptr;
			for (size_t i = 0; i < n_row; i++) {
				for (size_t j = data.rowptr[i]; j < data.rowptr[i + 1]; j++) {
					newcolptr_j[0] = i;
					s_copy(newcolptr_j + 1, nowcolptr_j, r - 1);
					newcolptr_j += r;
					nowcolptr_j += r - 1;
				}
			}
			s_free(data.colptr);
			data.colptr = newcolptr;

			data.rowptr = { 0, data.rowptr.back() };
			data.dims = prepend_num(data.dims, (size_t)1);
			data.rank++;
		}

		sparse_tensor(const sparse_mat<T, index_t>& mat, thread_pool* pool = nullptr) {
			*this = sparse_tensor<T, index_t, SPARSE_CSR>(mat, pool);
		}

		sparse_tensor& operator=(const sparse_mat<T, index_t>& mat) {
			*this = sparse_tensor<T, index_t, SPARSE_CSR>(mat);
			return *this;
		}

		void print_test() {
			for (size_t j = 0; j < data.rowptr[1]; j++) {
				for (size_t k = 0; k < data.rank - 1; k++)
					std::cout << (size_t)(data.colptr[j * (data.rank - 1) + k]) << " ";
				std::cout << " : " << data.valptr[j] << std::endl;
			}
		}
	};

	// some other functions

	template <typename T, typename index_t>
	inline T* sparse_vec_entry(const sparse_vec<T, index_t>& vec, const index_t index, const bool isbinary = true) {
		if (vec.nnz() == 0)
			return nullptr;
		index_t* ptr;
		if (isbinary)
			ptr = SparseRREF::binary_search(vec.indices, vec.indices + vec.nnz(), index);
		else
			ptr = std::find(vec.indices, vec.indices + vec.nnz(), index);
		if (ptr == vec.indices + vec.nnz())
			return nullptr;
		return vec.entries + (ptr - vec.indices);
	}

	template <typename T, typename index_t>
	inline T* sparse_mat_entry(sparse_mat<T, index_t>& mat, size_t r, index_t c, bool isbinary = true) {
		return sparse_vec_entry(mat[r], c, isbinary);
	}

	// join two sparse matrices
	template <typename T, typename index_t>
	sparse_mat<T, index_t> sparse_mat_join(const sparse_mat<T, index_t>& A, const sparse_mat<T, index_t>& B, thread_pool* pool = nullptr) {
		sparse_mat<T, index_t> res(A.nrow + B.nrow, std::max(A.ncol, B.ncol));

		if (pool == nullptr) {
			std::copy(A.rows.begin(), A.rows.end(), res.rows.begin());
			std::copy(B.rows.begin(), B.rows.end(), res.rows.begin() + A.nrow);
		}
		else {
			pool->detach_loop(0, A.nrow, [&](size_t i) {
				res[i] = A[i];
				});
			pool->wait();
			pool->detach_loop(0, B.nrow, [&](size_t i) {
				res[i + A.nrow] = B[i];
				});
			pool->wait();
		}

		return res;
	}

	template <typename T, typename index_t>
	sparse_mat<T, index_t> sparse_mat_join(sparse_mat<T, index_t>&& A, sparse_mat<T, index_t>&& B) {
		sparse_mat<T, index_t> res = std::move(A);
		res.append(std::move(B));
		res.ncol = std::max(A.ncol, B.ncol);
		return res;
	}

	// split a sparse matrix into two parts
	template <typename T, typename index_t>
	std::pair<sparse_mat<T, index_t>, sparse_mat<T, index_t>> sparse_mat_split(const sparse_mat<T, index_t>& mat, const size_t split_row, thread_pool* pool = nullptr) {
		if (split_row > mat.nrow)
			throw std::out_of_range("sparse_mat_split: split_row out of range");

		sparse_mat<T, index_t> A(split_row, mat.ncol);
		sparse_mat<T, index_t> B(mat.nrow - split_row, mat.ncol);

		if (pool == nullptr) {
			std::copy(mat.rows.begin(), mat.rows.begin() + split_row, A.rows.begin());
			std::copy(mat.rows.begin() + split_row, mat.rows.end(), B.rows.begin());
		}
		else {
			pool->detach_loop(0, split_row, [&](size_t i) {
				A[i] = mat[i];
				});
			pool->wait();
			pool->detach_loop(split_row, mat.nrow, [&](size_t i) {
				B[i - split_row] = mat[i];
				});
			pool->wait();
		}

		return { A, B };
	}

	template <typename T, typename index_t>
	std::pair<sparse_mat<T, index_t>, sparse_mat<T, index_t>> sparse_mat_split(sparse_mat<T, index_t>&& mat, const size_t split_row) {
		if (split_row > mat.nrow)
			throw std::out_of_range("sparse_mat_split: split_row out of range");

		sparse_mat<T, index_t> A(split_row, mat.ncol);
		sparse_mat<T, index_t> B(mat.nrow - split_row, mat.ncol);

		A.rows = std::vector<sparse_vec<T, index_t>>(std::make_move_iterator(mat.rows.begin()),
			std::make_move_iterator(mat.rows.begin() + split_row));
		B.rows = std::vector<sparse_vec<T, index_t>>(std::make_move_iterator(mat.rows.begin() + split_row),
			std::make_move_iterator(mat.rows.end()));

		mat.clear();

		return { A, B };
	}

	// a submatrix view of sparse_mat
	// it only contains a pointer to the original matrix and a list of row indices
	// and it does not own any memory
	// to save memory, if rows[0] > mat.nrow, then it is a full view
	template <typename T, typename index_t>
	struct sparse_mat_subview {
		sparse_mat<T, index_t>* mat_ptr = nullptr;
		std::vector<size_t> rows;

		sparse_mat_subview() = default;
		~sparse_mat_subview() = default;

		sparse_mat_subview(sparse_mat<T, index_t>& mat_) {
			mat_ptr = &mat_;
			rows = { mat_.nrow + 1 }; // full view
		}

		sparse_mat_subview(const sparse_mat<T, index_t>& mat_) {
			mat_ptr = const_cast<sparse_mat<T, index_t>*>(&mat_);
			rows = { mat_.nrow + 1 }; // full view
		}

		sparse_mat_subview(sparse_mat<T, index_t>& mat_, const std::vector<size_t>& rows_) {
			mat_ptr = &mat_;
			rows = rows_;
		}

		sparse_mat_subview(const sparse_mat_subview& l) { mat_ptr = l.mat_ptr; rows = l.rows; }
		sparse_mat_subview(sparse_mat_subview&& l) noexcept { mat_ptr = l.mat_ptr; rows = std::move(l.rows); }
		sparse_mat_subview& operator=(const sparse_mat_subview& l) {
			if (this == &l)
				return *this;
			mat_ptr = l.mat_ptr;
			rows = l.rows;
			return *this;
		}
		sparse_mat_subview& operator=(sparse_mat_subview&& l) noexcept {
			if (this == &l)
				return *this;
			mat_ptr = l.mat_ptr;
			rows = std::move(l.rows);
			return *this;
		}

		size_t nrow() const {
			if (mat_ptr) {
				if (rows[0] > mat_ptr->nrow) // full view
					return mat_ptr->nrow;
				return rows.size();
			}
			return 0;
		}

		size_t ncol() const {
			if (mat_ptr)
				return mat_ptr->ncol;
			return 0;
		}

		size_t nnz() const {
			if (mat_ptr) {
				size_t result = 0;
				if (rows[0] > mat_ptr->nrow) // full view
					return mat_ptr->nnz();
				for (auto r : rows)
					result += (*mat_ptr)[r].nnz();
				return result;
			}
			return 0;
		}

		bool is_full() const {
			if (mat_ptr) {
				return rows[0] > mat_ptr->nrow; // full view
			}
			return false;
		}

		// access the i-th row of the subview, it is dangerous because we do not check the index
		sparse_vec<T, index_t>& operator[](size_t i) {
			if (rows[0] > mat_ptr->nrow) // full view
				return (*mat_ptr)[i];
			return (*mat_ptr)[rows[i]];
		}
		const sparse_vec<T, index_t>& operator[](size_t i) const {
			if (rows[0] > mat_ptr->nrow) // full view
				return (*mat_ptr)[i];
			return (*mat_ptr)[rows[i]];
		}

		const size_t operator()(size_t i) const {
			if (rows[0] > mat_ptr->nrow) // full view
				return i;
			return rows[i];
		}

		sparse_mat<T, index_t>& get_mat() {
			return *mat_ptr;
		}
		const sparse_mat<T, index_t>& get_mat() const {
			return *mat_ptr;
		}

		void traverse(std::function<void(size_t)> func) {
			if (mat_ptr == nullptr)
				return;
			if (rows[0] > mat_ptr->nrow) { // full view
				for (size_t i = 0; i < mat_ptr->nrow; i++)
					func(i);
			}
			else {
				for (auto r : rows)
					func(r);
			}
		}
	};

}

#endif