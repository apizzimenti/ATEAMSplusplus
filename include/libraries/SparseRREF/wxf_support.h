/*
	Copyright (C) 2025 Zhenjie Li (Li, Zhenjie)

	This file is part of SparseRREF. The SparseRREF is free software:
	you can redistribute it and/or modify it under the terms of the MIT
	License.
*/

#ifndef WXF_SUPPORT_H
#define WXF_SUPPORT_H

#include "wxf_parser.h"
#include "sparse_type.h"

namespace SparseRREF {
	// TODO: reuse the code for sparse_mat_read_wxf and sparse_tensor_read_wxf

	template <typename T, typename index_t>
	sparse_mat<T, index_t> sparse_mat_read_wxf(const std::vector<WXF_PARSER::Token>& tokens, const field_t& F) {
		if (tokens.size() == 0)
			return sparse_mat<T, index_t>();

		if (tokens[0].type != WXF_PARSER::WXF_HEAD::func ||
			tokens[0].length != 4 ||
			tokens[1].type != WXF_PARSER::WXF_HEAD::symbol ||
			tokens[1].get_string_view() != "SparseArray" ||
			tokens[2].type != WXF_PARSER::WXF_HEAD::symbol ||
			tokens[2].get_string_view() != "Automatic") {
			std::cerr << "Error: sparse_mat_read: not a SparseArray" << std::endl;
			return sparse_mat<T, index_t>();
		}

#define GENERATE_COPY_ARR(TYPE, CTYPE, FUNC) \
		case TYPE: { \
			auto sp = token.get_arr_span<CTYPE>(); \
			for (const auto& v : sp) out.push_back(FUNC(v)); \
			break; }

#define TMP_IDENTITY_FUNC(x) (x)

#define GENERATE_ALL_ARR(FUNC2) \
		GENERATE_COPY_ARR(0, int8_t, FUNC2)  \
		GENERATE_COPY_ARR(1, int16_t, FUNC2) \
		GENERATE_COPY_ARR(2, int32_t, FUNC2) \
		GENERATE_COPY_ARR(3, int64_t, FUNC2) \
		GENERATE_COPY_ARR(16, uint8_t, FUNC2)  \
		GENERATE_COPY_ARR(17, uint16_t, FUNC2) \
		GENERATE_COPY_ARR(18, uint32_t, FUNC2) \
		GENERATE_COPY_ARR(19, uint64_t, FUNC2)

		auto copy_arr = [](const WXF_PARSER::Token& token, auto& out) {
			int num_type = token.dimensions[0];
			out.reserve(token.dimensions[1]);
			switch (num_type) {
				GENERATE_ALL_ARR(TMP_IDENTITY_FUNC);
			default:
				std::cerr << "Error: sparse_mat_read: read array fails" << std::endl;
				break;
			}
			};

		auto copy_modify_arr = [](const WXF_PARSER::Token& token, auto& out, auto&& func) {
			int num_type = token.dimensions[0];
			out.reserve(token.dimensions[1]);
			switch (num_type) {
				GENERATE_ALL_ARR(func);
			default:
				std::cerr << "Error: sparse_mat_read: read array fails" << std::endl;
				break;
			}
			};
#undef TMP_IDENTITY_FUNC
#undef GENERATE_COPY_ARR
#undef GENERATE_ALL_ARR

		// dims
		std::vector<size_t> dims;
		copy_arr(tokens[3], dims);

		// imp_val
		if (tokens[4].get_integer() != 0) {
			std::cerr << "Error: sparse_mat_read: the implicit value is not 0" << std::endl;
			return sparse_mat<T, index_t>();
		}

		// List of length 3
		// {1,{rowptr,colindex},vals}

		// check 
		if (tokens[5].type != WXF_PARSER::WXF_HEAD::func ||
			tokens[5].length != 3 ||
			tokens[6].type != WXF_PARSER::WXF_HEAD::symbol ||
			tokens[6].get_string_view() != "List" ||
			tokens[8].type != WXF_PARSER::WXF_HEAD::func ||
			tokens[8].length != 2 ||
			tokens[9].type != WXF_PARSER::WXF_HEAD::symbol ||
			tokens[9].get_string_view() != "List") {
			std::cerr << "Error: sparse_mat_read: wrong format in SparseArray" << std::endl;
			return sparse_mat<T, index_t>();
		}

		// rowptr is tokens[10]
		std::vector<size_t> rowptr;
		copy_arr(tokens[10], rowptr);

		// colindex is tokens[11]
		std::vector<index_t> colindex;
		if (WXF_PARSER::size_of_arr_num_type(tokens[11].dimensions[0]) > sizeof(index_t)) {
			std::cerr << "Error: sparse_mat_read: the type of index is not enough for colindex" << std::endl;
			return sparse_mat<T, index_t>();
		}
		copy_arr(tokens[11], colindex);

		std::vector<char> buffer;
		auto get_int_from_sv = [&buffer](std::string_view sv) -> int_t {
			buffer.clear();
			buffer.reserve(sv.size() + 1);

			buffer.insert(buffer.end(), sv.begin(), sv.end());
			buffer.push_back('\0');
			return int_t(buffer.data());
			};

		auto get_int_from_tv = [&get_int_from_sv](const WXF_PARSER::Token& node) -> int_t {
			switch (node.type) {
			case WXF_PARSER::WXF_HEAD::i8:
			case WXF_PARSER::WXF_HEAD::i16:
			case WXF_PARSER::WXF_HEAD::i32:
			case WXF_PARSER::WXF_HEAD::i64:
				return int_t(node.get_integer());
			case WXF_PARSER::WXF_HEAD::bigint:
				return get_int_from_sv(node.get_string_view());
			default:
				std::cerr << "not a integer" << std::endl;
				return int_t(0);
			}
			};

		// the other is vals
		std::vector<T> vals;
		if (tokens[12].type == WXF_PARSER::WXF_HEAD::array ||
			tokens[12].type == WXF_PARSER::WXF_HEAD::narray) {
			if constexpr (std::is_same_v<T, rat_t>)
				copy_arr(tokens[12], vals);
			else if constexpr (std::is_same_v<T, ulong>) {
				copy_modify_arr(tokens[12], vals, [&](auto v) {
					return nmod_set_si(v, F.mod);
					});
			}
		}
		else {
			// it is a list
			if (tokens[12].type != WXF_PARSER::WXF_HEAD::func ||
				tokens[13].type != WXF_PARSER::WXF_HEAD::symbol ||
				tokens[13].get_string_view() != "List") {
				std::cerr << "Error: sparse_mat_read: wrong format in SparseArray" << std::endl;
				return sparse_mat<T, index_t>();
			}

			auto nnz = tokens[12].length;
			vals.reserve(nnz);
			size_t pos = 14;
			while (pos < tokens.size()) {
				auto& token = tokens[pos];
				T val;
				switch (token.type) {
				case WXF_PARSER::WXF_HEAD::i8:
				case WXF_PARSER::WXF_HEAD::i16:
				case WXF_PARSER::WXF_HEAD::i32:
				case WXF_PARSER::WXF_HEAD::i64:
					if constexpr (std::is_same_v<T, rat_t>) {
						val = token.get_integer();
					}
					else if constexpr (std::is_same_v<T, ulong>) {
						val = int_t(token.get_integer()) % F.mod;
					}
					break;
				case WXF_PARSER::WXF_HEAD::bigint:
					if constexpr (std::is_same_v<T, rat_t>) {
						val = get_int_from_tv(token);
					}
					else if constexpr (std::is_same_v<T, ulong>) {
						val = get_int_from_sv(token.get_string_view()) % F.mod;
					}
					break;
				case WXF_PARSER::WXF_HEAD::func: {
					auto& ntoken = tokens[pos + 1];
					if (ntoken.get_string_view() == "Rational") {
						if (pos + 3 >= tokens.size()) {
							std::cerr << "Error: sparse_mat_read: wrong format in SparseArray" << std::endl;
							return sparse_mat<T, index_t>();
						}

						int_t n_1 = get_int_from_tv(tokens[pos + 2]);
						int_t d_1 = get_int_from_tv(tokens[pos + 3]);
						if constexpr (std::is_same_v<T, rat_t>) {
							val = rat_t(std::move(n_1), std::move(d_1), true);
						}
						else if constexpr (std::is_same_v<T, ulong>) {
							val = rat_t(std::move(n_1), std::move(d_1), true) % F.mod;
						}
						pos += 3;
					}
					else {
						std::cerr << "Error: sparse_mat_read: ";
						std::cerr << "not a SparseArray with rational / integer entries" << std::endl;
						return sparse_mat<T, index_t>();
					}
					break;
				}
				default:
					std::cerr << "Error: sparse_mat_read: ";
					std::cerr << "not a SparseArray with rational / integer entries" << std::endl;
					return sparse_mat<T, index_t>();
					break;
				}
				vals.push_back(val);
				pos++;
			}
		}

		sparse_mat<T, index_t> mat(dims[0], dims[1]);
		for (size_t i = 0; i < rowptr.size() - 1; i++) {
			mat[i].reserve(rowptr[i + 1] - rowptr[i]);
			for (auto j = rowptr[i]; j < rowptr[i + 1]; j++) {
				// mathematica is 1-indexed
				mat[i].push_back(colindex[j] - 1, vals[j]);
			}
		}

		return mat;
	}

	template <typename T, typename index_t>
	sparse_mat<T, index_t> sparse_mat_read_wxf(const std::filesystem::path file, const field_t& F) {
		auto fz = std::filesystem::file_size(file);

		WXF_PARSER::Parser wxf_parser;

		// if > 1GB, use mmap
		if (fz > 1ULL << 30) {
			MMapFile mm;
			std::string cc_str = std::filesystem::canonical(file).string();
			bool status = mmap_file(cc_str.c_str(), mm);

			if (!status) {
				// if mmap failed, read the file directly
				auto buffer = file_to_ustr(file);
				wxf_parser = WXF_PARSER::Parser(buffer);
			}

			wxf_parser = WXF_PARSER::Parser(mm.view);
		}

		auto buffer = file_to_ustr(file);
		wxf_parser = WXF_PARSER::Parser(buffer);

		wxf_parser.parse();

		return sparse_mat_read_wxf<T, index_t>(wxf_parser.tokens, F);
	}

	// SparseArray[Automatic,dims,imp_val = 0,{1,{rowptr,colindex},vals}]
	template <typename T, typename index_t>
	std::vector<uint8_t> sparse_mat_write_wxf(const sparse_mat<T, index_t>& mat, bool include_head = true) {
		using namespace WXF_PARSER;

		if (mat.nrow == 0)
			return std::vector<uint8_t>();

		std::string_view ff_template = "SparseArray[Automatic,#dims,0,{1,{#rowptr,#colindex},#vals}]";

		size_t rank = 2;
		size_t nnz = mat.nnz();
		std::vector<int64_t> dims_arr = { (int64_t)mat.nrow, (int64_t)mat.ncol };

		std::unordered_map<std::string, std::function<void(Encoder&)>> func_map;

		func_map["#dims"] = [&](Encoder& enc) {
			enc.push_packed_array({ 2 }, dims_arr);
			};

		func_map["#rowptr"] = [&](Encoder& enc) {
			std::vector<int64_t> rowptr(mat.nrow + 1);
			rowptr[0] = 0;
			for (size_t i = 0; i < mat.nrow; i++) {
				rowptr[i + 1] = rowptr[i] + mat[i].nnz();
			}
			enc.push_packed_array({ mat.nrow + 1 }, rowptr);
			};

		func_map["#colindex"] = [&](Encoder& enc) {
			auto mz = minimal_pos_signed_bits(mat.ncol);
			enc.push_array_info({ nnz, rank - 1 }, WXF_HEAD::array, mz);

#define APPEND_COLIND_DATA(TYPE) do {                                   \
		std::vector<TYPE> tmp((rank - 1) * nnz);                        \
		size_t idx = 0;                                                 \
		for (size_t i = 0; i < mat.nrow; i++) {                         \
			for (size_t j = 0; j < mat[i].nnz(); j++) {                 \
				tmp[idx] = mat[i](j) + 1;                               \
				idx++;                                                  \
			}                                                           \
		}                                                               \
		enc.push_ustr(tmp.data(), tmp.size());                          \
	} while (0)

			switch (mz) {
			case 0: APPEND_COLIND_DATA(int8_t); break;
			case 1: APPEND_COLIND_DATA(int16_t); break;
			case 2: APPEND_COLIND_DATA(int32_t); break;
			case 3: APPEND_COLIND_DATA(int64_t); break;
			case 4: std::cerr << "Error: sparse_tensor_write_wxf: too large dimension" << std::endl; break;
			}
#undef APPEND_COLIND_DATA
			};

		func_map["#vals"] = [&](Encoder& enc) {
			auto push_int = [&](const int_t& val) {
				if (val.fits_si())
					enc.push_integer(val.to_si());
				else
					enc.push_bigint(val.get_str());
				};

			if constexpr (std::is_same_v<T, rat_t>) {
				enc.push_function("List", nnz);
				for (size_t i = 0; i < mat.nrow; i++) {
					for (size_t j = 0; j < mat[i].nnz(); j++) {
						auto& rat_val = mat[i][j];
						if (rat_val.is_integer())
							push_int(rat_val.num());
						else {
							// func,2,symbol,8,"Rational"
							enc.push_ustr("f\x02s\x08Rational");
							push_int(rat_val.num());
							push_int(rat_val.den());
						}
					}
				}
			}
			else if constexpr (std::is_same_v<T, int_t>) {
				enc.push_function("List", nnz);
				for (size_t i = 0; i < mat.nrow; i++) {
					for (size_t j = 0; j < mat[i].nnz(); j++)
						push_int(mat[i][j]);
				}
			}
			else if constexpr (std::is_same_v<T, ulong>) {
				enc.push_array_info({ nnz }, WXF_HEAD::narray, 19);
				for (size_t i = 0; i < mat.nrow; i++) {
					enc.push_ustr(mat[i].entries, mat[i].nnz());
				}
			}
			};

		Encoder enc = fullform_to_wxf(ff_template, func_map, include_head);

		return enc.buffer;
	}

	template<typename T, typename IndexType, typename S>
	void sparse_tensor_write(S& st, const sparse_tensor<T, IndexType, SPARSE_COO>& tensor) {
		const auto& dims = tensor.dims();
		const size_t rank = dims.size();
		char num_buf[32];

		for (size_t i = 0; i < rank; ++i) {
			auto [ptr, ec] = std::to_chars(num_buf, num_buf + sizeof(num_buf), dims[i]);
			st.write(num_buf, ptr - num_buf);
			st.put(' ');
		}
		auto [ptr, ec] = std::to_chars(num_buf, num_buf + sizeof(num_buf), tensor.nnz());
		st.write(num_buf, ptr - num_buf);
		st.put('\n');

		std::vector<char> index_buf;
		index_buf.reserve(rank * 20 + 64);

		for (size_t i = 0; i < tensor.nnz(); ++i) {
			index_buf.clear();
			const auto& index = tensor.index(i);
			for (size_t j = 0; j < rank; ++j) {
				auto [ptr, ec] = std::to_chars(num_buf, num_buf + sizeof(num_buf), index[j] + 1);
				index_buf.insert(index_buf.end(), num_buf, ptr);
				index_buf.push_back(' ');
			}
			st.write(index_buf.data(), index_buf.size());
			st << tensor.val(i) << '\n';
		}
	}

	// SparseArray[Automatic,dims,imp_val = 0,{1,{rowptr,colindex},vals}]
	template <typename T, typename index_t>
	sparse_tensor<T, index_t, SPARSE_CSR> sparse_tensor_read_wxf(const std::vector<WXF_PARSER::Token>& tokens, const field_t& F, thread_pool* pool = nullptr, bool sort_ind = true) {
		using st = sparse_tensor<T, index_t, SPARSE_CSR>;

		if (tokens.size() == 0)
			return st();

		if (tokens[0].type != WXF_PARSER::WXF_HEAD::func ||
			tokens[0].length != 4 ||
			tokens[1].type != WXF_PARSER::WXF_HEAD::symbol ||
			tokens[1].get_string_view() != "SparseArray" ||
			tokens[2].type != WXF_PARSER::WXF_HEAD::symbol ||
			tokens[2].get_string_view() != "Automatic") {
			std::cerr << "Error: sparse_tensor_read: not a SparseArray" << std::endl;
			return st();
		}

#define GENERATE_COPY_ARR(TYPE, CTYPE, FUNC) \
		case TYPE: { \
			auto sp = token.get_arr_span<CTYPE>(); \
			for (const auto& v : sp) {*out = FUNC(v); out++;} \
			break; }

#define TMP_IDENTITY_FUNC(x) (x)

#define GENERATE_ALL_ARR(FUNC2) \
		GENERATE_COPY_ARR(0, int8_t, FUNC2)  \
		GENERATE_COPY_ARR(1, int16_t, FUNC2) \
		GENERATE_COPY_ARR(2, int32_t, FUNC2) \
		GENERATE_COPY_ARR(3, int64_t, FUNC2) \
		GENERATE_COPY_ARR(16, uint8_t, FUNC2)  \
		GENERATE_COPY_ARR(17, uint16_t, FUNC2) \
		GENERATE_COPY_ARR(18, uint32_t, FUNC2) \
		GENERATE_COPY_ARR(19, uint64_t, FUNC2)

		auto copy_arr = [](const WXF_PARSER::Token& token, auto out) {
			int num_type = token.dimensions[0];
			switch (num_type) {
				GENERATE_ALL_ARR(TMP_IDENTITY_FUNC);
			default:
				std::cerr << "Error: sparse_tensor_read: read array fails" << std::endl;
				break;
			}
			};

		auto copy_modify_arr = [](const WXF_PARSER::Token& token, auto out, auto&& func) {
			int num_type = token.dimensions[0];
			switch (num_type) {
				GENERATE_ALL_ARR(func);
			default:
				std::cerr << "Error: sparse_tensor_read: read array fails" << std::endl;
				break;
			}
			};
#undef TMP_IDENTITY_FUNC
#undef GENERATE_COPY_ARR
#undef GENERATE_ALL_ARR

		// dims
		std::vector<size_t> dims(tokens[3].dimensions[1]);
		copy_arr(tokens[3], dims.data());

		// imp_val
		if (tokens[4].get_integer() != 0) {
			std::cerr << "Error: sparse_tensor_read: the implicit value is not 0" << std::endl;
			return st();
		}

		// List of length 3
		// {1,{rowptr,colindex},vals}

		// check 
		if (tokens[5].type != WXF_PARSER::WXF_HEAD::func ||
			tokens[5].length != 3 ||
			tokens[6].type != WXF_PARSER::WXF_HEAD::symbol ||
			tokens[6].get_string_view() != "List" ||
			tokens[8].type != WXF_PARSER::WXF_HEAD::func ||
			tokens[8].length != 2 ||
			tokens[9].type != WXF_PARSER::WXF_HEAD::symbol ||
			tokens[9].get_string_view() != "List") {
			std::cerr << "Error: sparse_tensor_read: wrong format in SparseArray" << std::endl;
			return st();
		}

		// rowptr is tokens[10]
		std::vector<size_t> rowptr(tokens[10].dimensions[1]);
		copy_arr(tokens[10], rowptr.data());
		size_t nz = rowptr.back();

		st tensor(dims, nz);
		tensor.data.rowptr = std::move(rowptr);

		// colindex is tokens[11]
		if (WXF_PARSER::size_of_arr_num_type(tokens[11].dimensions[0]) > sizeof(index_t)) {
			std::cerr << "Error: sparse_tensor_read: the type of index is not enough for colindex" << std::endl;
			return st();
		}
		// mma is 1-based for colindex
		copy_modify_arr(tokens[11], tensor.data.colptr, [](auto v) { return v - 1; });

		std::vector<char> buffer;
		auto get_int_from_sv = [&buffer](std::string_view sv) -> int_t {
			buffer.clear();
			buffer.reserve(sv.size() + 1);

			buffer.insert(buffer.end(), sv.begin(), sv.end());
			buffer.push_back('\0');
			return int_t(buffer.data());
			};

		auto get_int_from_tv = [&get_int_from_sv](const WXF_PARSER::Token& node) -> int_t {
			switch (node.type) {
			case WXF_PARSER::WXF_HEAD::i8:
			case WXF_PARSER::WXF_HEAD::i16:
			case WXF_PARSER::WXF_HEAD::i32:
			case WXF_PARSER::WXF_HEAD::i64:
				return int_t(node.get_integer());
			case WXF_PARSER::WXF_HEAD::bigint:
				return get_int_from_sv(node.get_string_view());
			default:
				std::cerr << "not a integer" << std::endl;
				return int_t(0);
			}
			};

		// the other is vals
		if (tokens[12].type == WXF_PARSER::WXF_HEAD::array ||
			tokens[12].type == WXF_PARSER::WXF_HEAD::narray) {
			if constexpr (std::is_same_v<T, rat_t>)
				copy_arr(tokens[12], tensor.data.valptr);
			else if constexpr (std::is_same_v<T, ulong>) {
				copy_modify_arr(tokens[12], tensor.data.valptr, [&](auto v) {
					return nmod_set_si(v, F.mod);
					});
			}
		}
		else {
			// it is a list
			if (tokens[12].type != WXF_PARSER::WXF_HEAD::func ||
				tokens[13].type != WXF_PARSER::WXF_HEAD::symbol ||
				tokens[13].get_string_view() != "List") {
				std::cerr << "Error: sparse_tensor_read: wrong format in SparseArray" << std::endl;
				return st();
			}

			size_t pos = 14;
			T* vals = tensor.data.valptr;
			while (pos < tokens.size()) {
				auto& token = tokens[pos];
				T val;
				switch (token.type) {
				case WXF_PARSER::WXF_HEAD::i8:
				case WXF_PARSER::WXF_HEAD::i16:
				case WXF_PARSER::WXF_HEAD::i32:
				case WXF_PARSER::WXF_HEAD::i64:
					if constexpr (std::is_same_v<T, rat_t>) {
						val = token.get_integer();
					}
					else if constexpr (std::is_same_v<T, ulong>) {
						val = int_t(token.get_integer()) % F.mod;
					}
					break;
				case WXF_PARSER::WXF_HEAD::bigint:
					if constexpr (std::is_same_v<T, rat_t>) {
						val = get_int_from_tv(token);
					}
					else if constexpr (std::is_same_v<T, ulong>) {
						val = get_int_from_sv(token.get_string_view()) % F.mod;
					}
					break;
				case WXF_PARSER::WXF_HEAD::func: {
					auto ntoken = tokens[pos + 1];
					if (ntoken.get_string_view() == "Rational") {
						if (pos + 3 >= tokens.size()) {
							std::cerr << "Error: sparse_tensor_read: wrong format in SparseArray" << std::endl;
							return sparse_mat<T, index_t>();
						}

						int_t n_1 = get_int_from_tv(tokens[pos + 2]);
						int_t d_1 = get_int_from_tv(tokens[pos + 3]);
						if constexpr (std::is_same_v<T, rat_t>) {
							val = rat_t(std::move(n_1), std::move(d_1), true);
						}
						else if constexpr (std::is_same_v<T, ulong>) {
							val = rat_t(std::move(n_1), std::move(d_1), true) % F.mod;
						}
						pos += 3;
					}
					else {
						std::cerr << "Error: sparse_tensor_read: ";
						std::cerr << "not a SparseArray with rational / integer entries" << std::endl;
						return st();
					}
					break;
				}
				default:
					std::cerr << "Error: sparse_tensor_read: ";
					std::cerr << "not a SparseArray with rational / integer entries" << std::endl;
					return st();
					break;
				}
				*vals = val;
				vals++;
				pos++;
			}
		}
		
		if (sort_ind && !tensor.check_sorted())
			tensor.sort_indices(pool);

		return tensor;
	}

	template <typename T, typename index_t>
	auto sparse_tensor_read_wxf(const std::filesystem::path file, const field_t& F, thread_pool* pool = nullptr, bool sort_ind = true) {
		auto fz = std::filesystem::file_size(file);

		WXF_PARSER::Parser wxf_parser;

		// if > 1GB, use mmap
		if (fz > 1ULL << 30) {
			MMapFile mm;
			std::string cc_str = std::filesystem::canonical(file).string();
			bool status = mmap_file(cc_str.c_str(), mm);

			if (!status) {
				// if mmap failed, read the file directly
				auto buffer = file_to_ustr(file);
				wxf_parser = WXF_PARSER::Parser(buffer);
			}

			wxf_parser = WXF_PARSER::Parser(mm.view);
		}

		auto buffer = file_to_ustr(file);
		wxf_parser = WXF_PARSER::Parser(buffer);

		wxf_parser.parse();

		return sparse_tensor_read_wxf<T, index_t>(wxf_parser.tokens, F, pool, sort_ind);
	}

	// SparseArray[Automatic,dims,imp_val = 0,{1,{rowptr,colindex},vals}]
	template <typename T, typename index_t>
	std::vector<uint8_t> sparse_tensor_write_wxf(const sparse_tensor<T, index_t, SPARSE_CSR>& tensor, bool include_head = true) {
		using namespace WXF_PARSER;

		if (tensor.alloc() == 0)
			return std::vector<uint8_t>();

		std::string_view ff_template = "SparseArray[Automatic,#dims,0,{1,{#rowptr,#colindex},#vals}]";

		auto rank = tensor.rank();
		auto nnz = tensor.nnz();
		std::vector<int64_t> dims(tensor.dims().begin(), tensor.dims().end());

		std::unordered_map<std::string, std::function<void(Encoder&)>> func_map;

		func_map["#dims"] = [&](Encoder& enc) {
			enc.push_packed_array({ rank }, dims);
			};

		// important: we assume size_t fits into int64_t !!!!!!!!!
		func_map["#rowptr"] = [&](Encoder& enc) {
			const auto& rowptr = tensor.data.rowptr;
			enc.push_array_info({ rowptr.size() }, WXF_HEAD::array, 3);
			enc.push_ustr(rowptr.data(), rowptr.size());
			};


		func_map["#colindex"] = [&](Encoder& enc) {
			auto mz = minimal_pos_signed_bits(1 + *std::max_element(dims.begin() + 1, dims.end()));
			enc.push_array_info({ nnz, rank - 1 }, WXF_HEAD::array, mz);

#define APPEND_COLIND_DATA(TYPE) do {                                   \
        std::vector<TYPE> tmp((rank - 1) * nnz);                        \
        for (size_t i = 0; i < tmp.size(); i++)                         \
            tmp[i] = tensor.data.colptr[i] + 1;                         \
		enc.push_ustr(tmp.data(), tmp.size());                          \
    } while (0)

			switch (mz) {
			case 0: APPEND_COLIND_DATA(int8_t); break;
			case 1: APPEND_COLIND_DATA(int16_t); break;
			case 2: APPEND_COLIND_DATA(int32_t); break;
			case 3: APPEND_COLIND_DATA(int64_t); break;
			case 4: std::cerr << "Error: sparse_tensor_write_wxf: too large dimension" << std::endl; break;
			}
#undef APPEND_COLIND_DATA
			};

		func_map["#vals"] = [&](Encoder& enc) {
			auto push_int = [&](const int_t& val) {
				if (val.fits_si())
					enc.push_integer(val.to_si());
				else
					enc.push_bigint(val.get_str());
				};

			if constexpr (std::is_same_v<T, rat_t>) {
				enc.push_function("List", nnz);
				for (size_t i = 0; i < nnz; i++) {
					auto& rat_val = tensor.val(i);
					if (rat_val.is_integer())
						push_int(rat_val.num());
					else {
						// func,2,symbol,8,"Rational"
						enc.push_ustr("f\x02s\x08Rational");
						push_int(rat_val.num());
						push_int(rat_val.den());
					}
				}
			}
			else if constexpr (std::is_same_v<T, int_t>) {
				enc.push_function("List", nnz);
				for (size_t i = 0; i < nnz; i++) {
					push_int(tensor.val(i));
				}
			}
			else if constexpr (std::is_same_v<T, ulong>) {
				enc.push_array_info({ nnz }, WXF_HEAD::narray, 19);
				enc.push_ustr(tensor.data.valptr, tensor.data.valptr + nnz);
			}
			};

		Encoder enc = fullform_to_wxf(ff_template, func_map, include_head);

		return enc.buffer;
	}
}

#endif // WXF_SUPPORT_H
