/*
	Copyright (C) 2025 Zhenjie Li (Li, Zhenjie)

	You can redistribute it and/or modify it under the terms of the MIT
	License.
*/

/*
	WXF is a binary format for faithfully serializing Wolfram Language expressions
	in a form suitable for outside storage or interchange with other programs.
	WXF can readily be interpreted using low-level native types available in many
	programming languages, making it suitable as a format for reading and writing
	Wolfram Language expressions in other programming languages.

	The details of the WXF format are described in the Wolfram Language documentation:
	https://reference.wolfram.com/language/tutorial/WXFFormatDescription.html.en .

	We here intend to support import and export a SparseArray expression with rational
	entries, so some types are not supported, such as complex numbers.
	The full list of supported types is given below:

	done	byte value  type of part
	*		102			function
	*		67			int8_t
	*		106			int16_t
	*		105			int32_t
	*		76			int64_t
	*		114			machine reals
	*		83			string
	*		66			binary string
	*		115			symbol
	*		73			big integer
	*		82			big real
	-		193			packed array
	-		194			numeric array
	*		65			association
	*		58			delayed rule in association
	*		45			rule in association

	* is supported, - is partially supported
*/

#pragma once

#include <iostream>
#include <filesystem>
#include <fstream>
#include <functional>
#include <vector>
#include <string>
#include <string_view>
#include <cstdint>
#include <cstring>
#include <variant>

namespace WXF_PARSER {

	enum class WXF_HEAD {
		// function type
		func = 102,
		association = 65,
		delay_rule = 58,
		rule = 45,
		// string type
		symbol = 115,
		string = 83,
		binary_string = 66,
		bigint = 73,
		bigreal = 82,
		// number type
		i8 = 67,
		i16 = 106,
		i32 = 105,
		i64 = 76,
		f64 = 114,
		// array type
		array = 193,
		narray = 194
	};

	using NumberType = std::variant<int8_t, int16_t, int32_t, int64_t,
		uint8_t, uint16_t, uint32_t, uint64_t,
		float, double, bool>;

	NumberType select_type(int index) {
		switch (index) {
		case 0: return int8_t(0);
		case 1: return int16_t(0);
		case 2: return int32_t(0);
		case 3: return int64_t(0);
		case 16: return uint8_t(0);
		case 17: return uint16_t(0);
		case 18: return uint32_t(0);
		case 19: return uint64_t(0);
		case 34: return float(0);
		case 35: return double(0);
		default:
			std::cerr << "Unsupported type index" << std::endl;
			return bool(0);
		}
	}

	template <typename T>
	requires std::is_integral_v<T>&& std::is_signed_v<T>
	inline uint8_t minimal_signed_bits(T x) noexcept {
		if (x >= INT8_MIN && x <= INT8_MAX) return 0;
		if (x >= INT16_MIN && x <= INT16_MAX) return 1;
		if (x >= INT32_MIN && x <= INT32_MAX) return 2;
		return 3; // for int64_t
	}

	// for positive signed/unsigned integer
	template <typename T>
	inline uint8_t minimal_pos_signed_bits(T x) noexcept {
		if (x <= INT8_MAX) return 0;
		if (x <= INT16_MAX) return 1;
		if (x <= INT32_MAX) return 2;
		if (x <= INT64_MAX) return 3;
		return 4; // for uint64_t
	}

	template <typename T>
	requires std::is_integral_v<T>&& std::is_unsigned_v<T>
	inline uint8_t minimal_unsigned_bits(T x) noexcept {
		if (x <= UINT8_MAX) return 0;
		if (x <= UINT16_MAX) return 1;
		if (x <= UINT32_MAX) return 2;
		return 3; // for uint64_t
	}

	inline void serialize_varint(std::vector<uint8_t>& buffer, uint64_t val) {
		while (val > 0) {
			uint8_t byte = val & 127;
			val >>= 7;
			if (val > 0) byte |= 128; // set the continuation bit
			buffer.push_back(byte);
		}
	}

	template <typename T>
	void serialize_binary(std::vector<uint8_t>& buffer, const T& value) {
		const size_t old_size = buffer.size();
		buffer.resize(old_size + sizeof(T));
		std::memcpy(buffer.data() + old_size, &value, sizeof(T));
	}

	template <typename T>
	void serialize_binary(std::vector<uint8_t>& buffer, const T* valptr, const size_t len) {
		const size_t old_size = buffer.size();
		buffer.resize(old_size + sizeof(T) * len);
		std::memcpy(buffer.data() + old_size, valptr, sizeof(T) * len);
	}

	struct TOKEN {
		WXF_HEAD type;
		int rank = 0;
		union {
			// for number, string, symbol, bigint
			uint64_t length;
			// for array and narray, dimensions[0] is the type, dimensions[1] is the total flatten length
			// so the length is dimensions is rank + 2
			uint64_t* dimensions;
		};
		union { // data
			int64_t i;
			double d;
			int64_t* i_arr; // for array
			uint64_t* u_arr; // only for narray
			double* d_arr; // for array and narray, but not fully supported yet
			char* str;
		};

		TOKEN() : type(WXF_HEAD::i8), rank(0), length(0), i(0) {}

		uint64_t dim(int i) const {
			if (rank > 0)
				return dimensions[i + 2];
			return length;
		}

		void clear() {
			if (type == WXF_HEAD::symbol
				|| type == WXF_HEAD::bigint
				|| type == WXF_HEAD::bigreal
				|| type == WXF_HEAD::string
				|| type == WXF_HEAD::binary_string) {
				free(str);
			}
			else if (type == WXF_HEAD::array || type == WXF_HEAD::narray) {
				free(i_arr);
				free(dimensions);
			}
			// no need to clear i, length, rank, type, as they are just basic types
		}

		~TOKEN() { clear(); }

		// disable copy constructor and copy assignment operator
		TOKEN(const TOKEN&) = delete;
		TOKEN& operator=(const TOKEN&) = delete;

		// move constructor
		TOKEN(TOKEN&& other) noexcept : type(other.type), rank(other.rank), length(other.length), i(other.i) {
			if (type == WXF_HEAD::symbol
				|| type == WXF_HEAD::bigint
				|| type == WXF_HEAD::bigreal
				|| type == WXF_HEAD::string
				|| type == WXF_HEAD::binary_string) {
				str = other.str;
				other.str = nullptr;
			}
			else if (type == WXF_HEAD::array || type == WXF_HEAD::narray) {
				dimensions = other.dimensions;
				i_arr = other.i_arr; // for array, since it is union, we can use i_arr for narray
				other.dimensions = nullptr;
				other.i_arr = nullptr;
			}
		}

		// symbol/bigint/string/binary_string, length, str
		TOKEN(WXF_HEAD t, const std::string_view s) : type(t), rank(0) {
			length = s.size();
			str = (char*)malloc(length + 1);
			std::memcpy(str, s.data(), length);
			str[length] = '\0'; // null-terminate the string
		}

		TOKEN(WXF_HEAD t, const char* s, const size_t len) : type(t), rank(0) {
			length = len;
			str = (char*)malloc(length + 1);
			std::memcpy(str, s, length);
			str[length] = '\0'; // null-terminate the string
		}

		// machine number, val (length is given by the sizeof(val))
		TOKEN(WXF_HEAD t, int8_t v) : type(t), rank(0), length(1), i(v) {}
		TOKEN(WXF_HEAD t, int16_t v) : type(t), rank(0), length(2), i(v) {}
		TOKEN(WXF_HEAD t, int32_t v) : type(t), rank(0), length(3), i(v) {}
		TOKEN(WXF_HEAD t, int64_t v) : type(t), rank(0), length(4), i(v) {}
		TOKEN(WXF_HEAD t, float v) : type(t), rank(0), length(2), d(v) {}
		TOKEN(WXF_HEAD t, double v) : type(t), rank(0), length(4), d(v) {}

		// function, association, delay_rule, rule
		TOKEN(WXF_HEAD t, uint64_t len) : type(t), rank(0), length(len), i(0) {}

		// array/narray
		// with_arr: whether to allocate the i_arr/u_arr, default true
		//	         otherwise, the user may want to use the existing data pointer, 
		//           or only use it for dimensions
		TOKEN(WXF_HEAD t, const std::vector<size_t>& dims, int num_type, size_t len, bool with_arr = true) : type(t) {
			int r = dims.size();
			rank = r;
			dimensions = (uint64_t*)malloc((r + 2) * sizeof(uint64_t));
			dimensions[0] = num_type;
			dimensions[1] = len;
			for (auto i = 0; i < r; i++) {
				dimensions[i + 2] = dims[i];
			}
			if (with_arr)
				i_arr = (int64_t*)malloc(len * sizeof(int64_t));
			else
				i_arr = nullptr;
		}

		void to_ustr(std::vector<uint8_t>& res) const {
			const auto& token = *this;

			switch (token.type) {
			case WXF_HEAD::i8:
			case WXF_HEAD::i16:
			case WXF_HEAD::i32:
			case WXF_HEAD::i64: {
				auto num_type = minimal_signed_bits(token.i);
				switch (num_type) {
				case 0:
					res.push_back((uint8_t)WXF_HEAD::i8);
					serialize_binary(res, (int8_t)token.i);
					break;
				case 1:
					res.push_back((uint8_t)WXF_HEAD::i16);
					serialize_binary(res, (int16_t)token.i);
					break;
				case 2:
					res.push_back((uint8_t)WXF_HEAD::i32);
					serialize_binary(res, (int32_t)token.i);
					break;
				case 3:
					res.push_back((uint8_t)WXF_HEAD::i64);
					serialize_binary(res, token.i);
					break;
				default:
					break;
				}
				break;
			}
			case WXF_HEAD::f64:
				res.push_back((uint8_t)WXF_HEAD::f64);
				serialize_binary(res, token.d);
				break;
			case WXF_HEAD::func:
			case WXF_HEAD::association: {
				res.push_back((uint8_t)token.type);
				serialize_varint(res, token.length);
				break;
			}
			case WXF_HEAD::rule:
			case WXF_HEAD::delay_rule:
				res.push_back((uint8_t)token.type);
				break;
			case WXF_HEAD::symbol:
			case WXF_HEAD::bigint:
			case WXF_HEAD::bigreal:
			case WXF_HEAD::string:
			case WXF_HEAD::binary_string:
				res.push_back((uint8_t)token.type);
				serialize_varint(res, token.length);
				res.insert(res.end(), token.str, token.str + token.length);
				break;
			case WXF_HEAD::array:
			case WXF_HEAD::narray:
				res.push_back((uint8_t)token.type);
				res.push_back(token.dimensions[0]);
				serialize_varint(res, token.rank);
				for (auto i = 0; i < token.rank; i++) {
					serialize_varint(res, token.dimensions[i + 2]);
				}
				if (token.i_arr == nullptr)
					break;
				if (token.dimensions[0] == 3) {
					serialize_binary(res, token.i_arr, token.dimensions[1]);
				}
				else if (token.type == WXF_HEAD::narray && token.dimensions[0] == 19) {
					serialize_binary(res, token.u_arr, token.dimensions[1]);
				}
				else if (token.dimensions[0] == 35) {
					serialize_binary(res, token.d_arr, token.dimensions[1]);
				}
				else {
					std::cerr << "Unsupported number type in packed array or numeric array. " << std::endl;
				}
				break;
			default:
				break;
			}
		}

		std::vector<uint8_t> to_ustr() const {
			std::vector<uint8_t> res;
			to_ustr(res);
			return res;
		}

		template<typename T>
		void print(T& ss) const {
			auto& token = *this;
			switch (token.type) {
			case WXF_HEAD::i8:
				ss << "i8: " << token.i << std::endl;
				break;
			case WXF_HEAD::i16:
				ss << "i16: " << token.i << std::endl;
				break;
			case WXF_HEAD::i32:
				ss << "i32: " << token.i << std::endl;
				break;
			case WXF_HEAD::i64:
				ss << "i64: " << token.i << std::endl;
				break;
			case WXF_HEAD::f64:
				ss << "f64: " << token.d << std::endl;
				break;
			case WXF_HEAD::symbol:
				ss << "symbol: " << token.str << std::endl;
				break;
			case WXF_HEAD::bigint:
				ss << "bigint: " << token.str << std::endl;
				break;
			case WXF_HEAD::bigreal:
				ss << "bigreal: " << token.str << std::endl;
				break;
			case WXF_HEAD::string:
				ss << "string: " << token.str << std::endl;
				break;
			case WXF_HEAD::binary_string:
				ss << "binary_string: " << token.str << std::endl;
				break;
			case WXF_HEAD::func:
				ss << "func: " << token.length << " vars" << std::endl;
				break;
			case WXF_HEAD::association:
				ss << "association: " << token.length << " rules" << std::endl;
				break;
			case WXF_HEAD::delay_rule:
				ss << "delay_rule: " << token.length << std::endl;
				break;
			case WXF_HEAD::rule:
				ss << "rule: " << token.length << std::endl;
				break;
			case WXF_HEAD::array: {
				ss << "array: rank = " << token.rank << ", dimensions = ";
				size_t all_len = token.dimensions[1];
				for (int i = 0; i < token.rank; i++) {
					ss << token.dimensions[i + 2] << " ";
				}
				ss << std::endl;

				auto num_type = token.dimensions[0];
				ss << "data: ";
				if (token.i_arr == nullptr)
					break;
				if (num_type < 4) {
					for (size_t i = 0; i < all_len; i++) {
						ss << token.i_arr[i] << " ";
					}
				}
				else if (num_type >= 34 && num_type <= 35) {
					for (size_t i = 0; i < all_len; i++) {
						ss << token.d_arr[i] << " ";
					}
				}
				else {
					std::cerr << "Unknown type" << std::endl;
				}
				ss << std::endl;
				break;
			}
			case WXF_HEAD::narray: {
				ss << "narray: rank = " << token.rank << ", dimensions = ";
				for (int i = 0; i < token.rank; i++) {
					ss << token.dimensions[i + 2] << " ";
				}
				ss << std::endl;

				size_t num_type = token.dimensions[0];
				size_t all_len = token.dimensions[1];

				ss << "data: ";
				if (token.i_arr == nullptr) 
					break;
				if (num_type >= 16 && num_type < 20) {
					for (size_t i = 0; i < all_len; i++)
						ss << token.u_arr[i] << " ";
				}
				else if (num_type < 4) {
					for (size_t i = 0; i < all_len; i++)
						ss << token.i_arr[i] << " ";
				}
				else if (num_type >= 34 && num_type <= 35) {
					for (size_t i = 0; i < all_len; i++)
						ss << token.d_arr[i] << " ";
				}
				else {
					std::cerr << "Unknown type" << std::endl;
				}
				ss << std::endl;
				break;
			}
			default:
				std::cerr << "Unknown type" << std::endl;
			}
		}

		void print() const {
			print(std::cout);
		}
	};

	struct Parser {
		const uint8_t* buffer; // the buffer to read
		size_t pos = 0;
		size_t size = 0; // the size of the buffer
		int err = 0; // 0 is ok, otherwise error
		std::vector<TOKEN> tokens; // the tokens read from the buffer

		Parser(const uint8_t* buf, const size_t len) : buffer(buf), pos(0), size(len), err(0) {}
		Parser(const std::vector<uint8_t>& buf) : buffer(buf.data()), pos(0), size(buf.size()), err(0) {}

		// we suppose that the length does not exceed 2^64 - 1 .. 
		uint64_t ReadVarint() {
			size_t count = 0;
			uint64_t result = 0;
			auto ptr = buffer + pos;

			while (pos < size && count < 8) {
				result |= (uint64_t)((*ptr) & 127) << (7 * count);
				count++; pos++;
				if (!((*ptr) & 128))
					break;
				ptr++;
			}

			return result;
		}

		void parseExpr() {
			// check the file head
			if (pos == 0) {
				if (size < 2 || buffer[0] != 56 || buffer[1] != 58) {
					std::cerr << "Invalid WXF file" << std::endl;
					err = 1;
					return;
				}
				pos = 2;
			}

			while (pos < size) {
				WXF_HEAD type = (WXF_HEAD)(buffer[pos]); pos++;

				if (pos == size)
					break;

				switch (type) {
				case WXF_HEAD::i8: {
					int8_t val;
					std::memcpy(&val, buffer + pos, sizeof(int8_t));
					pos += sizeof(int8_t) / sizeof(uint8_t);
					tokens.emplace_back(type, val);
					break;
				}
				case WXF_HEAD::i16: {
					int16_t val;
					std::memcpy(&val, buffer + pos, sizeof(int16_t));
					pos += sizeof(int16_t) / sizeof(uint8_t);
					tokens.emplace_back(type, val);
					break;
				}
				case WXF_HEAD::i32: {
					int32_t val;
					std::memcpy(&val, buffer + pos, sizeof(int32_t));
					pos += sizeof(int32_t) / sizeof(uint8_t);
					tokens.emplace_back(type, val);
					break;
				}
				case WXF_HEAD::i64: {
					int64_t val;
					std::memcpy(&val, buffer + pos, sizeof(int64_t));
					pos += sizeof(int64_t) / sizeof(uint8_t);
					tokens.emplace_back(type, val);
					break;
				}
				case WXF_HEAD::f64: {
					double val;
					std::memcpy(&val, buffer + pos, sizeof(double));
					pos += sizeof(double) / sizeof(uint8_t);
					tokens.emplace_back(type, val);
					break;
				}
				case WXF_HEAD::symbol:
				case WXF_HEAD::bigint:
				case WXF_HEAD::bigreal:
				case WXF_HEAD::string:
				case WXF_HEAD::binary_string: {
					auto length = ReadVarint();
					tokens.emplace_back(type, (const char*)(buffer + pos), length);
					pos += length;
					break;
				}
				case WXF_HEAD::func:
				case WXF_HEAD::association:
					tokens.emplace_back(type, uint64_t(ReadVarint()));
					break;
				case WXF_HEAD::delay_rule:
				case WXF_HEAD::rule:
					tokens.emplace_back(type, uint64_t(2));
					break;
				case WXF_HEAD::array:
				case WXF_HEAD::narray: {
					auto num_type = ReadVarint();
					if (num_type > 50) {
						std::cerr << "Unsupported type: " << num_type << std::endl;
						err = 2;
						break;
					}
					auto r = ReadVarint();
					std::vector<size_t> dims(r);
					size_t all_len = 1;
					for (size_t i = 0; i < r; i++) {
						dims[i] = ReadVarint();
						all_len *= dims[i];
					}
					tokens.emplace_back(type, dims, num_type, all_len);
					makeArray(type, tokens.back());
					break;
				}
				default:
					std::cerr << "Unknown head type: " << (int)type << " pos: " << pos << std::endl;
					err = 3;
					break;
				}
			}
			err = 0;
		}

		// numeric/packed array, rank, dimensions, data
		// for the num_type
		// 0 is int8_t      1 is int16_t
		// 2 is int32_t     3 is int64_t
		// 16 is uint8_t    17 is uint16_t ; only for numeric array
		// 18 is uint32_t   19 is uint64_t ; only for numeric array
		// 34 float         35 double
		// 51 complex float 52 complex double
		// we only support int8_t, int16_t, int32_t, int64_t, float, double

		void makeArray(WXF_HEAD type, TOKEN& token) {
			auto num_type = token.dimensions[0];
			auto all_len = token.dimensions[1];
			size_t size_of_type;
			if (num_type >= 16 && num_type < 20) {
				size_of_type = (size_t)1 << (num_type - 16);
			}
			else if (num_type >= 34 && num_type <= 35) {
				size_of_type = (size_t)1 << (num_type - 32);
			}
			else {
				size_of_type = (size_t)1 << num_type;
			}

			auto load_array = [&](auto&& assign) {
				std::visit([&](auto&& x) {
					using T = std::decay_t<decltype(x)>;
					T val;
					for (size_t i = 0; i < all_len; i++) {
						std::memcpy(&val, buffer + pos, size_of_type);
						assign(i, val);
						pos += size_of_type;
					}
					}, select_type(num_type));
				};

			if (type == WXF_HEAD::narray && num_type >= 16 && num_type < 20) {
				load_array([&](size_t i, auto val) { token.u_arr[i] = static_cast<uint64_t>(val); });
			}
			else if (num_type < 10) {
				load_array([&](size_t i, auto val) { token.i_arr[i] = static_cast<int64_t>(val); });
			}
			else if (num_type >= 34 && num_type <= 35) {
				load_array([&](size_t i, auto val) { token.d_arr[i] = static_cast<double>(val); });
			}
		}
	};

	struct ExprNode {
		size_t index; // the index of the token in the tokens vector
		size_t size; // the size of the children
		std::unique_ptr<ExprNode[]> children; // the children of the node
		WXF_HEAD type;

		ExprNode() : index(0), size(0), children(nullptr), type(WXF_HEAD::i8) {} // default constructor

		ExprNode(size_t idx, size_t sz, WXF_HEAD t) : index(idx), size(sz), type(t) {
			constexpr size_t MAX_ALLOC = std::numeric_limits<int64_t>::max();
			if (size > MAX_ALLOC) {
				throw std::bad_alloc();
			} else if (size > 0) {
				children = std::make_unique<ExprNode[]>(size);
			}
		}

		ExprNode(const ExprNode&) = delete; // disable copy constructor
		ExprNode& operator=(const ExprNode&) = delete; // disable copy assignment operator
		// move constructor
		ExprNode(ExprNode&& other) noexcept : index(other.index), size(other.size),
			children(std::move(other.children)), type(other.type) {
			other.index = 0;
			other.size = 0;
			other.children = nullptr;
			other.type = WXF_HEAD::i8;
		}


		// move assignment operator
		ExprNode& operator=(ExprNode&& other) noexcept {
			if (this != &other) {
				index = other.index;
				size = other.size;
				children = std::move(other.children);
				type = other.type;
				other.size = 0;
				other.index = 0;
				other.children = nullptr;
				other.type = WXF_HEAD::i8;
			}
			return *this;
		}

		// destructor
		void clear() {
			index = 0;
			size = 0;
			type = WXF_HEAD::i8;
			if (children) {
				children.reset();
			}
		}

		~ExprNode() {
			clear();
		}

		bool has_children() const {
			return size > 0;
		}

		const ExprNode& operator[] (size_t i) const {
			return children[i];
		}

		ExprNode& operator[] (size_t i) {
			return children[i];
		}
	};

	void node_to_ustr(std::vector<uint8_t>& res, const std::vector<TOKEN>& tokens, const ExprNode& node) {
		auto& token = tokens[node.index];

		switch (node.type) {
		case WXF_HEAD::func:
		case WXF_HEAD::association: 
			res.push_back((uint8_t)node.type);
			serialize_varint(res, node.size);
			res.push_back((uint8_t)token.type);
			serialize_varint(res, token.length);
			res.insert(res.end(), token.str, token.str + token.length);
			for (size_t i = 0; i < node.size; i++) {
				node_to_ustr(res, tokens, node.children[i]);
			}
			break;
		case WXF_HEAD::rule:
		case WXF_HEAD::delay_rule:
			res.push_back((uint8_t)node.type);
			for (size_t i = 0; i < node.size; i++) {
				node_to_ustr(res, tokens, node.children[i]);
			}
			break;
		default:
			token.to_ustr(res);
			break;
		}
	}

	struct ExprTree {
		std::vector<TOKEN> tokens;
		ExprNode root;

		ExprTree() {} // default constructor
		ExprTree(Parser parser, size_t index, size_t size, WXF_HEAD type) : root(index, size, type) {
			tokens = std::move(parser.tokens);
		}

		ExprTree(const ExprTree&) = delete; // disable copy constructor
		ExprTree& operator=(const ExprTree&) = delete; // disable copy assignment operator

		// move constructor
		ExprTree(ExprTree&& other) noexcept : tokens(std::move(other.tokens)), root(std::move(other.root)) {
			other.root.size = 0;
			other.root.index = 0;
			other.root.children = nullptr;
		}

		// move assignment operator
		ExprTree& operator=(ExprTree&& other) noexcept {
			if (this != &other) {
				tokens = std::move(other.tokens);
				root = std::move(other.root);
				other.root.size = 0;
				other.root.index = 0;
				other.root.children = nullptr;
			}
			return *this;
		}

		const TOKEN& operator[](const ExprNode& node) const {
			return tokens[node.index];
		}

		// it is super slow...
		// TODO: optimize this function
		std::vector<uint8_t> to_ustr(bool include_head = true) const {
			std::vector<uint8_t> res;
			if (include_head) {
				res.push_back(56); // WXF head
				res.push_back(58); // WXF head
			}
			node_to_ustr(res, tokens, root);
			return res;
		}

		//void plot() const {
		//	printPrettyTree(root);
		//}
	};

	ExprTree MakeExprTree(Parser& parser) {
		ExprTree tree;
		if (parser.err != 0)
			return tree;

		tree.tokens = std::move(parser.tokens);
		auto total_len = tree.tokens.size();
		auto& tokens = tree.tokens;

		std::vector<ExprNode*> expr_stack; // the stack to store the current father nodes
		std::vector<size_t> node_stack; // the vector to store the node index

		std::function<void(void)> move_to_next_node = [&]() {
			if (node_stack.empty())
				return;

			node_stack.back()++; // move to the next node
			if (node_stack.back() >= expr_stack.back()->size) {
				expr_stack.pop_back(); // pop the current node
				node_stack.pop_back(); // pop the current node index
				move_to_next_node();
			}
			};

		// first we need to find the root node
		size_t pos = 0;
		auto& token = tokens[pos];
		if (token.type == WXF_HEAD::func) {
			// i + 1 is the head of the function (a symbol)
			tree.root = ExprNode(pos + 1, token.length, token.type);
			pos += 2; // skip the head
		}
		else if (token.type == WXF_HEAD::association) {
			// association does not have a head
			tree.root = ExprNode(pos + 1, token.length, token.type);
			pos += 1;
		}
		else {
			// if the token is not a function type, only one token is allowed
			tree.root = ExprNode(pos, 0, token.type);
			return tree;
		}

		expr_stack.push_back(&(tree.root));
		node_stack.push_back(0);

		// now we need to parse the expression
		for (; pos < total_len; pos++) {
			auto& token = tokens[pos];
			if (token.type == WXF_HEAD::func || token.type == WXF_HEAD::association) {
				// if the token is a function type, we need to create a new node
				auto node_pos = node_stack.back();
				auto parent = expr_stack.back();
				auto& node = parent->children[node_pos];
				if (token.type == WXF_HEAD::func) {
					node = ExprNode(pos + 1, token.length, token.type);
					pos++; // skip the head
				}
				else
					node = ExprNode(pos, token.length, token.type);
				expr_stack.push_back(&(node)); // push the new node to the stack
				node_stack.push_back(0); // push the new node index to the stack
			}
			else if (token.type == WXF_HEAD::delay_rule || token.type == WXF_HEAD::rule) {
				// if the token is a rule type, we need to create a new node
				auto node_pos = node_stack.back();
				auto parent = expr_stack.back();
				auto& node = parent->children[node_pos];
				node = ExprNode(pos, 2, token.type);
				expr_stack.push_back(&(node)); // push the new node to the stack
				node_stack.push_back(0); // push the new node index to the stack
			}
			else {
				// if the token is not a function type, we need to move to the next node
				auto node_pos = node_stack.back();
				auto parent = expr_stack.back();
				auto& node = parent->children[node_pos];
				node = ExprNode(pos, 0, token.type);

				move_to_next_node();
			}
		}

		if (!node_stack.empty()) {
			std::cerr << "Error: not all nodes are parsed" << std::endl;
			for (auto& node : expr_stack) {
				node->clear();
			}
		}

		return tree;
	}

	ExprTree MakeExprTree(const uint8_t* str, const size_t len) {
		Parser parser(str, len);
		parser.parseExpr();
		return MakeExprTree(parser);
	}

	ExprTree MakeExprTree(const std::vector<uint8_t>& str) {
		Parser parser(str);
		parser.parseExpr();
		return MakeExprTree(parser);
	}

	ExprTree MakeExprTree(const std::string_view str) {
		Parser parser(reinterpret_cast<const uint8_t*>(str.data()), str.size());
		parser.parseExpr();
		return MakeExprTree(parser);
	}

	ExprTree MakeExprTree(const std::filesystem::path filename) {
		if (!std::filesystem::exists(filename)) {
			std::cerr << "Error: File does not exist!" << std::endl;
			return ExprTree();
		}
		ExprTree expr_tree;
		std::ifstream file(filename, std::ios::binary | std::ios::ate);
		std::streamsize size = file.tellg();
		file.seekg(0, std::ios::beg);

		std::vector<uint8_t> buffer(size);
		if (!file.read((char*)buffer.data(), size)) {
			std::cerr << "Failed to read file!" << std::endl;
			return ExprTree();
		}

		return MakeExprTree(buffer);
	}

	// debug only: convert the expression tree to DOT format of Graphviz
	template <typename SS>
	void toDotFormat(const ExprTree& tree, SS& oss) {
		oss << "digraph ExprTree {\n";
		oss << "  node [shape=box];\n";

		auto& tokens = tree.tokens;

		int nodeId = 0;
		std::function<void(const ExprNode&)> traverse = [&](const ExprNode& node) {
			int currentId = nodeId++;

			oss << "  n" << currentId << " [label=\"";
			switch (tokens[node.index].type) {
			case WXF_HEAD::symbol:
				oss << tokens[node.index].str;
				break;
			case WXF_HEAD::i8:
			case WXF_HEAD::i16:
			case WXF_HEAD::i32:
			case WXF_HEAD::i64:
				oss << tokens[node.index].i;
				break;
			case WXF_HEAD::f64:
				oss << tokens[node.index].d;
				break;
			case WXF_HEAD::bigint:
				oss << "bigint " << node.index;
				break;
			case WXF_HEAD::bigreal:
				oss << "bigreal " << node.index;
				break;
			case WXF_HEAD::string:
				oss << "string " << node.index;
				break;
			case WXF_HEAD::binary_string:
				oss << "binary_string " << node.index;
				break;
			case WXF_HEAD::array:
				oss << "array " << node.index;
				break;
			case WXF_HEAD::narray:
				oss << "narray " << node.index;
				break;
			default:
				oss << node.index;
				break;
			}
			oss << "\"];\n";

			for (size_t i = 0; i < node.size; ++i) {
				int childId = nodeId;
				traverse(node.children[i]);
				oss << "  n" << currentId << " -> n" << childId << ";\n";
			}
			};

		traverse(tree.root);
		oss << "}\n";
	}

	/*
		an example to test:
			SparseArray[{{1, 1} -> 1/3.0, {1, 23133} ->
			N[Pi, 100] + I N[E, 100], {44, 2} -> -(4/
			 33333333333333444333333335), {_, _} -> 0}]

		FullForm:
			SparseArray[Automatic,List[44,23133],0,
			List[1,List[List[0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
			2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3],
			List[List[1],List[23133],List[2]]],
			List[0.3333333333333333`,
			Complex[3.1415926535897932384626433832795028841971693993751058209
			749445923078164062862089986280348253421170679821480865191976`100.,
			2.7182818284590452353602874713526624977572470936999595749669676277
			240766303535475945713821785251664274274663919320031`100.],
			Rational[-4,33333333333333444333333335]]]]

		std::vector<uint8_t> test{ 56, 58, 102, 4, 115, 11, 83, 112, 97, 114, 115, 101, 65, 114, 114, \
								97, 121, 115, 9, 65, 117, 116, 111, 109, 97, 116, 105, 99, 193, 1, 1, \
								2, 44, 0, 93, 90, 67, 0, 102, 3, 115, 4, 76, 105, 115, 116, 67, 1, \
								102, 2, 115, 4, 76, 105, 115, 116, 193, 0, 1, 45, 0, 2, 2, 2, 2, 2, \
								2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, \
								2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 193, 1, 2, 3, 1, 1, \
								0, 93, 90, 2, 0, 102, 3, 115, 4, 76, 105, 115, 116, 114, 85, 85, 85, \
								85, 85, 85, 213, 63, 102, 2, 115, 7, 67, 111, 109, 112, 108, 101, \
								120, 82, 122, 51, 46, 49, 52, 49, 53, 57, 50, 54, 53, 51, 53, 56, 57, \
								55, 57, 51, 50, 51, 56, 52, 54, 50, 54, 52, 51, 51, 56, 51, 50, 55, \
								57, 53, 48, 50, 56, 56, 52, 49, 57, 55, 49, 54, 57, 51, 57, 57, 51, \
								55, 53, 49, 48, 53, 56, 50, 48, 57, 55, 52, 57, 52, 52, 53, 57, 50, \
								51, 48, 55, 56, 49, 54, 52, 48, 54, 50, 56, 54, 50, 48, 56, 57, 57, \
								56, 54, 50, 56, 48, 51, 52, 56, 50, 53, 51, 52, 50, 49, 49, 55, 48, \
								54, 55, 57, 56, 50, 49, 52, 56, 48, 56, 54, 53, 49, 57, 49, 57, 55, \
								54, 96, 49, 48, 48, 46, 82, 122, 50, 46, 55, 49, 56, 50, 56, 49, 56, \
								50, 56, 52, 53, 57, 48, 52, 53, 50, 51, 53, 51, 54, 48, 50, 56, 55, \
								52, 55, 49, 51, 53, 50, 54, 54, 50, 52, 57, 55, 55, 53, 55, 50, 52, \
								55, 48, 57, 51, 54, 57, 57, 57, 53, 57, 53, 55, 52, 57, 54, 54, 57, \
								54, 55, 54, 50, 55, 55, 50, 52, 48, 55, 54, 54, 51, 48, 51, 53, 51, \
								53, 52, 55, 53, 57, 52, 53, 55, 49, 51, 56, 50, 49, 55, 56, 53, 50, \
								53, 49, 54, 54, 52, 50, 55, 52, 50, 55, 52, 54, 54, 51, 57, 49, 57, \
								51, 50, 48, 48, 51, 49, 96, 49, 48, 48, 46, 102, 2, 115, 8, 82, 97, \
								116, 105, 111, 110, 97, 108, 67, 252, 73, 26, 51, 51, 51, 51, 51, 51, \
								51, 51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 51, 51, 51, 51, 51, 51, \
								51, 51, 53 };

		auto tree = MakeExprTree(test);
		for (const auto& token : tree.tokens) {
			token.print();
		};

		result:
			func: 4 vars
			symbol: SparseArray
			symbol: Automatic
			array: rank = 1, dimensions = 2
			data: 44 23133
			i8: 0
			func: 3 vars
			symbol: List
			i8: 1
			func: 2 vars
			symbol: List
			array: rank = 1, dimensions = 45
			data: 0 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3
			array: rank = 2, dimensions = 3 1
			data: 1 23133 2
			func: 3 vars
			symbol: List
			f64: 0.333333
			func: 2 vars
			symbol: Complex
			bigreal: 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865191976`100.
			bigreal: 2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320031`100.
			func: 2 vars
			symbol: Rational
			i8: -4
			bigint: 33333333333333444333333335
	*/

} // namespace WXF_PARSER