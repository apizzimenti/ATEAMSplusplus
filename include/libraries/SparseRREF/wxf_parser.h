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
	*		193			packed array
	*		194			numeric array
	*		65			association
	*		58			delayed rule in association
	*		45			rule in association
*/

#pragma once

#include <complex>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <functional>
#include <vector>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>

namespace WXF_PARSER {

	using complex_float_t = std::complex<float>;
	using complex_double_t = std::complex<double>;

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

	inline size_t size_of_head_num_type(const WXF_HEAD head) {
		switch (head) {
		case WXF_HEAD::i8:
			return sizeof(int8_t);
		case WXF_HEAD::i16:
			return sizeof(int16_t);
		case WXF_HEAD::i32:
			return sizeof(int32_t);
		case WXF_HEAD::i64:
			return sizeof(int64_t);
		case WXF_HEAD::f64:
			return sizeof(double);
		default:
			return 0;
		}
	}

	// array: head(numeric/packed array), num_type, rank, dimensions, data
	// for the num_type
	// 0 is int8_t      1 is int16_t
	// 2 is int32_t     3 is int64_t
	// 16 is uint8_t    17 is uint16_t ; only for numeric array
	// 18 is uint32_t   19 is uint64_t ; only for numeric array
	// 34 float         35 double
	// 51 complex float 52 complex double

	// only the last 3 bits are used to indicate the size
	constexpr size_t size_of_arr_num_type(const int num_type) {
		return size_t(1) << (num_type & 0b111);
	}

	template <typename T>
		requires std::is_integral_v<T>&& std::is_signed_v<T>
	constexpr uint8_t minimal_signed_bits(T x) noexcept {
		if (x >= INT8_MIN && x <= INT8_MAX) return 0;
		if (x >= INT16_MIN && x <= INT16_MAX) return 1;
		if (x >= INT32_MIN && x <= INT32_MAX) return 2;
		return 3; // for int64_t
	}

	// for positive signed/unsigned integer
	template <typename T>
	constexpr uint8_t minimal_pos_signed_bits(T x) noexcept {
		if (x <= INT8_MAX) return 0;
		if (x <= INT16_MAX) return 1;
		if (x <= INT32_MAX) return 2;
		if (x <= INT64_MAX) return 3;
		return 4; // for uint64_t
	}

	template <typename T>
		requires std::is_integral_v<T>&& std::is_unsigned_v<T>
	constexpr uint8_t minimal_unsigned_bits(T x) noexcept {
		if (x <= UINT8_MAX) return 0;
		if (x <= UINT16_MAX) return 1;
		if (x <= UINT32_MAX) return 2;
		return 3; // for uint64_t
	}

	inline void serialize_varint(std::vector<uint8_t>& buffer, uint64_t val) {
		uint8_t temp[10];
		size_t i = 0;

		do {
			temp[i] = val & 0x7F;
			val >>= 7;
			if (val != 0) temp[i] |= 0x80;
			++i;
		} while (val != 0);

		buffer.insert(buffer.end(), temp, temp + i);
	}

	template <typename T>
	void serialize_binary(std::vector<uint8_t>& buffer, const T& value) {
		const size_t old_size = buffer.size();
		buffer.resize(old_size + sizeof(T));
		std::memcpy(buffer.data() + old_size, &value, sizeof(T));
	}

	template <typename T>
	void serialize_binary(std::vector<uint8_t>& buffer, const T* valptr, const size_t len) {
		buffer.insert(buffer.end(), (uint8_t*)valptr, (uint8_t*)(valptr + len));
	}

	struct Encoder {
		std::vector<uint8_t> buffer;

		Encoder() = default;
		~Encoder() = default;
		Encoder(const Encoder&) = default;
		Encoder& operator=(const Encoder&) = default;
		Encoder(Encoder&&) = default;
		Encoder& operator=(Encoder&&) = default;

		void clear() { buffer.clear(); }

		// move from existing buffer
		Encoder(std::vector<uint8_t>&& buf) : buffer(std::move(buf)) {}

		// push ustr directly
		Encoder& push_ustr(const std::vector<uint8_t>& str) { buffer.insert(buffer.end(), str.begin(), str.end()); return *this; }
		Encoder& push_ustr(const std::string_view str) {
			buffer.insert(buffer.end(), (uint8_t*)str.data(), (uint8_t*)(str.data() + str.size())); return *this;
		}
		template<typename T>
		Encoder& push_ustr(const T* str_ptr, const size_t len) {
			buffer.insert(buffer.end(), (uint8_t*)str_ptr, (uint8_t*)(str_ptr + len)); return *this;
		}
		template<typename T>
		Encoder& push_ustr(T* start, T* end) {
			buffer.insert(buffer.end(), (uint8_t*)start, (uint8_t*)end); return *this;
		}

		Encoder& push_integer(const int64_t val) {
			auto num_type = minimal_signed_bits(val);
			switch (num_type) {
			case 0:
				buffer.push_back((uint8_t)WXF_HEAD::i8);
				serialize_binary(buffer, (int8_t)val);
				break;
			case 1:
				buffer.push_back((uint8_t)WXF_HEAD::i16);
				serialize_binary(buffer, (int16_t)val);
				break;
			case 2:
				buffer.push_back((uint8_t)WXF_HEAD::i32);
				serialize_binary(buffer, (int32_t)val);
				break;
			case 3:
				buffer.push_back((uint8_t)WXF_HEAD::i64);
				serialize_binary(buffer, val);
				break;
			default:
				break;
			}
			return *this;
		}

		Encoder& push_real(const double val) {
			buffer.push_back((uint8_t)WXF_HEAD::f64);
			serialize_binary(buffer, val);
			return *this;
		}

		// push string type struct: string/symbol/bigint/bigreal, default type is string
		Encoder& push_string(const std::string_view str, const WXF_HEAD type = WXF_HEAD::string) {
			buffer.push_back((uint8_t)type);
			serialize_varint(buffer, str.size());
			return push_ustr(str);
		}

		Encoder& push_symbol(const std::string_view sym) { return push_string(sym, WXF_HEAD::symbol); }
		Encoder& push_bigint(const std::string_view bigint_str) { return push_string(bigint_str, WXF_HEAD::bigint); }
		Encoder& push_bigreal(const std::string_view bigreal_str) { return push_string(bigreal_str, WXF_HEAD::bigreal); }
		Encoder& push_binary_string(const std::string_view bin_str) { return push_string(bin_str, WXF_HEAD::binary_string); }

		Encoder& push_function(const std::string_view head, const size_t num_vars) {
			buffer.push_back((uint8_t)WXF_HEAD::func);
			serialize_varint(buffer, num_vars);
			return push_string(head, WXF_HEAD::symbol);
		}

		Encoder& push_association(const size_t num_rules) {
			buffer.push_back((uint8_t)WXF_HEAD::association);
			serialize_varint(buffer, num_rules);
			return *this;
		}

		Encoder& push_rule() { buffer.push_back((uint8_t)WXF_HEAD::rule); return *this; }
		Encoder& push_delay_rule() { buffer.push_back((uint8_t)WXF_HEAD::delay_rule); return *this; }

		// return the total length of the array
		size_t push_array_info(const std::vector<size_t>& dimension_array, WXF_HEAD type, uint8_t num_type) {
			size_t all_len = 1;
			// [array_type, num_type, rank, dimensions...]
			buffer.push_back((uint8_t)type);
			buffer.push_back(num_type);
			serialize_varint(buffer, dimension_array.size());
			for (auto dim : dimension_array) {
				serialize_varint(buffer, dim);
				all_len *= dim;
			}
			return all_len;
		}

		template<typename T>
		Encoder& push_array(const std::vector<size_t>& dimension_array, const std::span<T> data, WXF_HEAD type, uint8_t num_type) {
			// backup current size
			size_t old_size = buffer.size();

			// [array_type, num_type, rank, dimensions..., data...]
			auto all_len = push_array_info(dimension_array, type, num_type);

			if (all_len != data.size()) {
				std::cerr << "Encoder::push_array: Data size does not match the dimension array." << std::endl;
				// restore buffer
				buffer.resize(old_size);
				return *this;
			}

			// push data
			return push_ustr(data.data(), data.size());
		}

		template<typename T>
			requires std::is_integral_v<T>&& std::is_signed_v<T>
		Encoder& push_packed_array(const std::vector<size_t>& dimension_array, const std::span<const T> data) {
			int num_type = minimal_signed_bits(std::numeric_limits<T>::max());
			return push_array(dimension_array, data, WXF_HEAD::array, num_type);
		}

		Encoder& push_packed_array(const std::vector<size_t>& dimension_array, const std::span<const float> data) {
			return push_array(dimension_array, data, WXF_HEAD::array, 34);
		}

		Encoder& push_packed_array(const std::vector<size_t>& dimension_array, const std::span<const double> data) {
			return push_array(dimension_array, data, WXF_HEAD::array, 35);
		}

		Encoder& push_packed_array(const std::vector<size_t>& dimension_array, const std::span<const complex_float_t> data) {
			return push_array(dimension_array, data, WXF_HEAD::array, 51);
		}

		Encoder& push_packed_array(const std::vector<size_t>& dimension_array, const std::span<const complex_double_t> data) {
			return push_array(dimension_array, data, WXF_HEAD::array, 52);
		}

		template<typename T>
		Encoder& push_packed_array(const std::vector<size_t>& dimension_array, const std::vector<T>& data) {
			return push_packed_array(dimension_array, std::span<const T>(data));
		}

		template<typename T>
			requires std::is_integral_v<T>
		Encoder& push_numeric_array(const std::vector<size_t>& dimension_array, const std::span<const T> data) {
			bool is_signed = std::is_signed_v<T>;
			int num_type;

			if (is_signed)
				num_type = minimal_signed_bits(std::numeric_limits<T>::max());
			else
				num_type = 16 + minimal_unsigned_bits(std::numeric_limits<T>::max());
			return push_array(dimension_array, data, WXF_HEAD::narray, num_type);
		}

		Encoder& push_numeric_array(const std::vector<size_t>& dimension_array, const std::span<const float> data) {
			return push_array(dimension_array, data, WXF_HEAD::narray, 34);
		}

		Encoder& push_numeric_array(const std::vector<size_t>& dimension_array, const std::span<const double> data) {
			return push_array(dimension_array, data, WXF_HEAD::narray, 35);
		}

		Encoder& push_numeric_array(const std::vector<size_t>& dimension_array, const std::span<const complex_float_t> data) {
			return push_array(dimension_array, data, WXF_HEAD::narray, 51);
		}

		Encoder& push_numeric_array(const std::vector<size_t>& dimension_array, const std::span<const complex_double_t> data) {
			return push_array(dimension_array, data, WXF_HEAD::narray, 52);
		}

		template<typename T>
		Encoder& push_numeric_array(const std::vector<size_t>& dimension_array, const std::vector<T>& data) {
			return push_numeric_array(dimension_array, std::span<const T>(data));
		}
	};

	struct Token {
		WXF_HEAD type;
		int rank = 0;
		union {
			// for number, string, symbol, bigint
			size_t length;
			// for array and narray, dimensions[0] is the type, dimensions[1] is the total flatten length
			// so the length is dimensions is rank + 2
			size_t* dimensions;
		};
		const uint8_t* data; // pointer to the data in the original buffer

		Token() : type(WXF_HEAD::i8), rank(0), length(0), data(nullptr) {}
		Token(const WXF_HEAD t, const size_t len, const uint8_t* d) : type(t), rank(0), length(len), data(d) {}
		Token(const WXF_HEAD t, const std::vector<size_t>& dims, const int num_type, const size_t len, const uint8_t* d) : type(t), data(d) {
			int r = dims.size();
			rank = r;
			dimensions = (size_t*)malloc((r + 2) * sizeof(size_t));
			dimensions[0] = num_type;
			dimensions[1] = len;
			for (auto i = 0; i < r; i++) {
				dimensions[i + 2] = dims[i];
			}
		}

		size_t dim(size_t i) const {
			if (rank > 0)
				return dimensions[i + 2];
			return length;
		}

		template<typename T> T* get_ptr() const { return (T*)data; }

		~Token() {
			if (type == WXF_HEAD::array || type == WXF_HEAD::narray)
				free(dimensions);
		}

		Token(const Token& other) : type(other.type), rank(other.rank), length(other.length), data(other.data) {
			if (type == WXF_HEAD::array || type == WXF_HEAD::narray) {
				dimensions = (size_t*)malloc((rank + 2) * sizeof(size_t));
				std::memcpy(dimensions, other.dimensions, (rank + 2) * sizeof(size_t));
			}
		}

		Token(Token&& other) noexcept : type(other.type), rank(other.rank), length(other.length), data(other.data) {
			if (type == WXF_HEAD::array || type == WXF_HEAD::narray) {
				dimensions = other.dimensions;
				other.dimensions = nullptr;
			}
		}

		Token& operator=(const Token& other) {
			if (this != &other) {
				type = other.type;
				rank = other.rank;
				length = other.length;
				data = other.data;
				if (type == WXF_HEAD::array || type == WXF_HEAD::narray) {
					dimensions = (size_t*)malloc((rank + 2) * sizeof(size_t));
					std::memcpy(dimensions, other.dimensions, (rank + 2) * sizeof(size_t));
				}
			}
			return *this;
		}

		Token& operator=(Token&& other) noexcept {
			if (this != &other) {
				type = other.type;
				rank = other.rank;
				length = other.length;
				data = other.data;
				if (type == WXF_HEAD::array || type == WXF_HEAD::narray) {
					dimensions = other.dimensions;
					other.dimensions = nullptr;
				}
			}
			return *this;
		}

		int64_t get_integer() const {
			if (type == WXF_HEAD::i8)
				return *(int8_t*)data;
			else if (type == WXF_HEAD::i16)
				return *(int16_t*)data;
			else if (type == WXF_HEAD::i32)
				return *(int32_t*)data;
			else if (type == WXF_HEAD::i64)
				return *(int64_t*)data;
			else
				return 0;
		}

		double get_real() const {
			if (type == WXF_HEAD::f64)
				return *(double*)data;
			else
				return 0;
		}

		std::string_view get_string_view() const {
			if (type == WXF_HEAD::symbol
				|| type == WXF_HEAD::bigint
				|| type == WXF_HEAD::bigreal
				|| type == WXF_HEAD::string
				|| type == WXF_HEAD::binary_string) {
				return std::string_view((const char*)data, length);
			}
			else
				return std::string_view();
		}

		template<typename T>
		std::span<const T> get_arr_span() const {
			if (type != WXF_HEAD::array && type != WXF_HEAD::narray)
				return std::span<const T>();
			return std::span<const T>((T*)data, dimensions[1]);
		}

		// debug only, print the token info
		template<typename T>
		void print(T& ss) const {
			auto& token = *this;
			switch (token.type) {
			case WXF_HEAD::i8:
				ss << "i8: " << token.get_integer() << std::endl;
				break;
			case WXF_HEAD::i16:
				ss << "i16: " << token.get_integer() << std::endl;
				break;
			case WXF_HEAD::i32:
				ss << "i32: " << token.get_integer() << std::endl;
				break;
			case WXF_HEAD::i64:
				ss << "i64: " << token.get_integer() << std::endl;
				break;
			case WXF_HEAD::f64:
				ss << "f64: " << token.get_real() << std::endl;
				break;
			case WXF_HEAD::symbol:
				ss << "symbol: " << token.get_string_view() << std::endl;
				break;
			case WXF_HEAD::bigint:
				ss << "bigint: " << token.get_string_view() << std::endl;
				break;
			case WXF_HEAD::bigreal:
				ss << "bigreal: " << token.get_string_view() << std::endl;
				break;
			case WXF_HEAD::string:
				ss << "string: " << token.get_string_view() << std::endl;
				break;
			case WXF_HEAD::binary_string:
				ss << "binary_string:" << token.get_string_view() << std::endl;
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
				if (token.data == nullptr)
					break;

#define PRINT_ARRAY_CASE(TYPE1, TYPE2) do {                      \
				for (size_t i = 0; i < all_len; i++) {           \
					ss << (TYPE2)(get_ptr<TYPE1>()[i]) << " ";   \
				}} while (0)

				switch (num_type) {
				case 0: PRINT_ARRAY_CASE(int8_t, int64_t); break;
				case 1: PRINT_ARRAY_CASE(int16_t, int64_t); break;
				case 2: PRINT_ARRAY_CASE(int32_t, int64_t); break;
				case 3: PRINT_ARRAY_CASE(int64_t, int64_t); break;
				case 34: PRINT_ARRAY_CASE(float, double); break;
				case 35: PRINT_ARRAY_CASE(double, double); break;
				case 51: PRINT_ARRAY_CASE(complex_float_t, complex_double_t); break;
				case 52: PRINT_ARRAY_CASE(complex_double_t, complex_double_t); break;
				default: std::cerr << "Unknown number type: " << num_type << std::endl; break;
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
				if (token.data == nullptr)
					break;
				switch (num_type) {
				case 0: PRINT_ARRAY_CASE(int8_t, int64_t); break;
				case 1: PRINT_ARRAY_CASE(int16_t, int64_t); break;
				case 2: PRINT_ARRAY_CASE(int32_t, int64_t); break;
				case 3: PRINT_ARRAY_CASE(int64_t, int64_t); break;
				case 16: PRINT_ARRAY_CASE(uint8_t, uint64_t); break;
				case 17: PRINT_ARRAY_CASE(uint16_t, uint64_t); break;
				case 18: PRINT_ARRAY_CASE(uint32_t, uint64_t); break;
				case 19: PRINT_ARRAY_CASE(uint64_t, uint64_t); break;
				case 34: PRINT_ARRAY_CASE(float, double); break;
				case 35: PRINT_ARRAY_CASE(double, double); break;
				case 51: PRINT_ARRAY_CASE(complex_float_t, complex_double_t); break;
				case 52: PRINT_ARRAY_CASE(complex_double_t, complex_double_t); break;
				default: std::cerr << "Unknown number type: " << num_type << std::endl; break;
				}
#undef PRINT_ARRAY_CASE
				ss << std::endl;
				break;
			}
			default:
				std::cerr << "Unknown type: " << (int)token.type << std::endl;
				break;
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
		std::vector<Token> tokens;

		Parser(const uint8_t* buf, const size_t len) : buffer(buf), pos(0), size(len), err(0) {}
		Parser(const std::vector<uint8_t>& buf) : buffer(buf.data()), pos(0), size(buf.size()), err(0) {}
		Parser(const std::string_view buf) : buffer((const uint8_t*)buf.data()), pos(0), size(buf.size()), err(0) {}

		// default special member functions
		Parser() = default;
		~Parser() = default;
		Parser(const Parser&) = default;
		Parser& operator=(const Parser&) = default;
		Parser(Parser&&) noexcept = default;
		Parser& operator=(Parser&&) noexcept = default;

		inline uint64_t read_varint() {
			const uint8_t* ptr = buffer + pos;
			const uint8_t* end = buffer + size;
			uint64_t result = 0;
			uint8_t b;

			if (ptr >= end) return 0;

			b = *ptr++; result = uint64_t(b & 0x7F);         if (!(b & 0x80) || ptr >= end) goto done;
			b = *ptr++; result |= uint64_t(b & 0x7F) << 7;   if (!(b & 0x80) || ptr >= end) goto done;
			b = *ptr++; result |= uint64_t(b & 0x7F) << 14;  if (!(b & 0x80) || ptr >= end) goto done;
			b = *ptr++; result |= uint64_t(b & 0x7F) << 21;  if (!(b & 0x80) || ptr >= end) goto done;
			b = *ptr++; result |= uint64_t(b & 0x7F) << 28;  if (!(b & 0x80) || ptr >= end) goto done;
			b = *ptr++; result |= uint64_t(b & 0x7F) << 35;  if (!(b & 0x80) || ptr >= end) goto done;
			b = *ptr++; result |= uint64_t(b & 0x7F) << 42;  if (!(b & 0x80) || ptr >= end) goto done;
			b = *ptr++; result |= uint64_t(b & 0x7F) << 49;  if (!(b & 0x80) || ptr >= end) goto done;
			b = *ptr++; result |= uint64_t(b & 0x7F) << 56;  if (!(b & 0x80) || ptr >= end) goto done;
			b = *ptr++; result |= uint64_t(b & 0x7F) << 63;
		done:
			pos = ptr - buffer;
			return result;
		}

		void parse() {
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
				case WXF_HEAD::i8:
				case WXF_HEAD::i16:
				case WXF_HEAD::i32:
				case WXF_HEAD::i64:
				case WXF_HEAD::f64: {
					auto length = size_of_head_num_type(type);
					tokens.emplace_back(type, length, buffer + pos);
					pos += length;
					break;
				}
				case WXF_HEAD::symbol:
				case WXF_HEAD::bigint:
				case WXF_HEAD::bigreal:
				case WXF_HEAD::string:
				case WXF_HEAD::binary_string: {
					auto length = read_varint();
					tokens.emplace_back(type, length, buffer + pos);
					pos += length;
					break;
				}
				case WXF_HEAD::func:
				case WXF_HEAD::association: {
					auto length = read_varint();
					tokens.emplace_back(type, length, buffer + pos);
					break;
				}
				case WXF_HEAD::delay_rule:
				case WXF_HEAD::rule:
					tokens.emplace_back(type, size_t(2), buffer + pos);
					break;
				case WXF_HEAD::array:
				case WXF_HEAD::narray: {
					int num_type = read_varint();
					auto r = read_varint();
					std::vector<size_t> dims(r);
					size_t all_len = 1;
					for (size_t i = 0; i < r; i++) {
						dims[i] = read_varint();
						all_len *= dims[i];
					}
					tokens.emplace_back(type, dims, num_type, all_len, buffer + pos);
					pos += all_len * size_of_arr_num_type(num_type);
					break;
				}
				default:
					std::cerr << "Unknown head type: " << (int)type << " pos: " << pos << std::endl;
					err = 2;
					break;
				}
			}
			err = 0;
		}
	};

	struct expr_node {
		size_t index; // the index of the token in the tokens vector
		std::vector<expr_node> children;
		WXF_HEAD type;

		expr_node() : index(0), children(), type(WXF_HEAD::i8) {} // default constructor
		~expr_node() = default;

		expr_node(const expr_node&) = default; // copy constructor
		expr_node& operator=(const expr_node&) = default; // copy assignment operator
		expr_node(expr_node&& other) = default; // move constructor
		expr_node& operator=(expr_node&& other) = default; // move assignment operator

		expr_node(size_t idx, size_t sz, WXF_HEAD t) : index(idx), type(t) {
			if (sz > 0) 
				children.resize(sz);
		}

		size_t size() const { return children.size(); }
		bool has_children() const { return children.size() > 0; }
		const expr_node& operator[] (size_t i) const { return children[i]; }
		expr_node& operator[] (size_t i) { return children[i]; }
	};

	struct expr_tree {
		std::vector<Token> tokens;
		expr_node root;

		expr_tree() {} // default constructor
		expr_tree(Parser parser, size_t index, size_t size, WXF_HEAD type) : root(index, size, type) {
			tokens = std::move(parser.tokens);
		}

		expr_tree(const expr_tree&) = default; //  copy constructor
		expr_tree& operator=(const expr_tree&) = default; // copy assignment operator
		expr_tree(expr_tree&&) noexcept = default; // move constructor
		expr_tree& operator=(expr_tree&&) noexcept = default; // move assignment operator
		~expr_tree() = default;

		const Token& operator[](const expr_node& node) const {
			return tokens[node.index];
		}

		void print(std::ostream& ss, const expr_node& node, const int level = 0) const {
			for (int i = 0; i < level; i++)
				ss << "  ";
			ss << "Node type: " << (int)node.type << ", index: " << node.index << ", size: " << node.size() << std::endl;
			for (size_t i = 0; i < node.size(); i++) {
				print(ss, node.children[i], level + 1);
			}
		}

		void print(std::ostream& ss) const {
			print(ss, root, 0);
		}

		void print() const {
			print(std::cout);
		}
	};

	inline expr_tree make_expr_tree(Parser& parser) {
		expr_tree tree;
		if (parser.err != 0)
			return tree;

		tree.tokens = std::move(parser.tokens);

		auto total_len = tree.tokens.size();
		auto& tokens = tree.tokens;

		std::vector<expr_node*> expr_stack; // the stack to store the current father nodes
		std::vector<size_t> node_stack; // the vector to store the node index

		std::function<void(void)> move_to_next_node = [&]() {
			if (node_stack.empty())
				return;

			node_stack.back()++; // move to the next node
			if (node_stack.back() >= expr_stack.back()->size()) {
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
			tree.root = expr_node(pos + 1, token.length, token.type);
			pos += 2; // skip the head
		}
		else if (token.type == WXF_HEAD::association) {
			// association does not have a head
			tree.root = expr_node(pos + 1, token.length, token.type);
			pos += 1;
		}
		else {
			// if the token is not a function type, only one token is allowed
			tree.root = expr_node(pos, 0, token.type);
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
					node = expr_node(pos + 1, token.length, token.type);
					pos++; // skip the head
				}
				else
					node = expr_node(pos, token.length, token.type);
				expr_stack.push_back(&(node)); // push the new node to the stack
				node_stack.push_back(0); // push the new node index to the stack
			}
			else if (token.type == WXF_HEAD::delay_rule || token.type == WXF_HEAD::rule) {
				// if the token is a rule type, we need to create a new node
				auto node_pos = node_stack.back();
				auto parent = expr_stack.back();
				auto& node = parent->children[node_pos];
				node = expr_node(pos, 2, token.type);
				expr_stack.push_back(&(node)); // push the new node to the stack
				node_stack.push_back(0); // push the new node index to the stack
			}
			else {
				// if the token is not a function type, we need to move to the next node
				auto node_pos = node_stack.back();
				auto parent = expr_stack.back();
				auto& node = parent->children[node_pos];
				node = expr_node(pos, 0, token.type);

				move_to_next_node();
			}
		}

		if (!node_stack.empty()) {
			std::cerr << "Error: not all nodes are parsed" << std::endl;
		}

		return tree;
	}

	inline expr_tree make_expr_tree(const uint8_t* str, const size_t len) {
		Parser parser(str, len);
		parser.parse();
		return make_expr_tree(parser);
	}

	inline expr_tree make_expr_tree(const std::vector<uint8_t>& str) {
		Parser parser(str);
		parser.parse();
		return make_expr_tree(parser);
	}

	inline expr_tree make_expr_tree(const std::string_view str) {
		Parser parser(str);
		parser.parse();
		return make_expr_tree(parser);
	}

} // namespace WXF_PARSER

/***********************************************************************************/

// a simple FullForm parser, we only support [0-9,a-z,A-Z,$] in symbol names
// it is used to write a template engine to generate WXF files
// and we add a special expression type starting with # for subexpression labels
// e.g. #x for a integer x, e.g. #1, #2, ...
// and then we can replace these labels with actual expressions when generating WXF files

// !! do not use for high percision numbers or very large integers, it only supports standard C++ number formats
// !! it is not efficient and robust enough now, make your template simple
// !! and handle complicated sub-expressions manually by using push_** in Encoder

namespace WXF_PARSER::FullForm {
	enum class atom_type {
		Integer,
		Real,
		String,
		Symbol,
		Expression,
		Null
	};

	struct atom_expression {
		atom_type type_;
		std::string value_;

		atom_expression(atom_type type, const std::string_view value)
			: type_(type), value_(value) {
		}

		const atom_type get_type() const { return type_; }
		const std::string& get_value() const { return value_; }

		std::string to_FullForm() const {
			switch (type_) {
			case atom_type::String:
				return "\"" + value_ + "\"";
			default:
				return value_;
			}
		}
	};

	struct lexer {
		std::string input_;
		size_t position_;
		size_t length_;

		char current_char() const {
			return position_ < length_ ? input_[position_] : '\0';
		}

		void advance() {
			if (position_ < length_) position_++;
		}

		void skip_ws() {
			while (position_ < length_ && std::isspace(static_cast<unsigned char>(current_char()))) {
				advance();
			}
		}

		enum token_type {
			IDENTIFIER,    // 
			EXPRESSION,    // start with #, use as a label of subexpression, it is not standard mathematica 
			INTEGER,       // 
			REAL,          // 
			STRING,        // 
			LBRACKET,      // [
			RBRACKET,      // ]
			COMMA,         // ,
			END            // 
		};

		struct token {
			token_type type;
			std::string value;
			size_t position;
		};

		lexer(const std::string_view input)
			: input_(input), position_(0), length_(input.length()) {
		}

		token nextToken() {
			skip_ws();

			if (position_ >= length_) {
				return { END, "", position_ };
			}

			size_t startPos = position_;
			char ch = current_char();

			// $ 
			if (std::isalpha(ch) || ch == '$') {
				std::string value;
				while (position_ < length_ &&
					(std::isalnum(static_cast<unsigned char>(current_char())) || current_char() == '$')) {
					value += current_char();
					advance();
				}
				return { IDENTIFIER, value, startPos };
			}

			// expression starting with #
			if (ch == '#') {
				std::string value;
				value += ch;
				advance();
				while (position_ < length_ &&
					(std::isalnum(static_cast<unsigned char>(current_char())) || current_char() == '$')) {
					value += current_char();
					advance();
				}
				return { EXPRESSION, value, startPos };
			}

			// number
			if (std::isdigit(static_cast<unsigned char>(ch)) || ch == '.' || ch == '-') {
				std::string value;
				bool hasDot = false;
				bool hasE = false;

				// minus sign
				if (ch == '-') {
					value += ch;
					advance();
					if (position_ >= length_ || !std::isdigit(static_cast<unsigned char>(current_char()))) {
						// single minus sign, return as IDENTIFIER
						return { IDENTIFIER, value, startPos };
					}
					ch = current_char();
				}

				// integer part
				while (position_ < length_ && std::isdigit(static_cast<unsigned char>(current_char()))) {
					value += current_char();
					advance();
				}

				// decimal part
				if (position_ < length_ && current_char() == '.') {
					hasDot = true;
					value += '.';
					advance();
					while (position_ < length_ && std::isdigit(static_cast<unsigned char>(current_char()))) {
						value += current_char();
						advance();
					}
				}

				// scientific notation
				if (position_ < length_ && (current_char() == 'e' || current_char() == 'E')) {
					hasE = true;
					value += current_char();
					advance();

					if (position_ < length_ && (current_char() == '+' || current_char() == '-')) {
						value += current_char();
						advance();
					}

					if (position_ >= length_ || !std::isdigit(static_cast<unsigned char>(current_char()))) {
						std::cerr << "Invalid scientific notation at position " << position_ << std::endl;
						return { END, "", startPos };
					}

					while (position_ < length_ && std::isdigit(static_cast<unsigned char>(current_char()))) {
						value += current_char();
						advance();
					}
				}

				if (hasDot || hasE) {
					return { REAL, value, startPos };
				}
				else {
					return { INTEGER, value, startPos };
				}
			}

			// string
			if (ch == '"') {
				std::string value;
				advance();

				while (position_ < length_ && current_char() != '"') {
					if (current_char() == '\\' && position_ + 1 < length_) {
						advance();
						char next = current_char();
						switch (next) {
						case 'n': value += '\n'; break;
						case 't': value += '\t'; break;
						case 'r': value += '\r'; break;
						case '"': value += '"'; break;
						case '\\': value += '\\'; break;
						default: value += '\\'; value += next; break;
						}
					}
					else {
						value += current_char();
					}
					advance();
				}

				if (position_ >= length_ || current_char() != '"') {
					std::cerr << "Unterminated string at position " << startPos << std::endl;
					return { END, "", startPos };
				}
				advance();  // skip closing "

				return { STRING, value, startPos };
			}

			// single-character tokens
			switch (ch) {
			case '[':
				advance();
				return { LBRACKET, "[", startPos };
			case ']':
				advance();
				return { RBRACKET, "]", startPos };
			case ',':
				advance();
				return { COMMA, ",", startPos };
			}

			std::cerr << "Unknown character '" << ch << "' at position " << position_ << std::endl;
			return { END, "", position_ };
		}
	};

	struct expression {
		atom_expression head_;
		std::vector<expression> args_;

		expression(const atom_expression& head) : head_(head), args_() {}
		expression(const atom_expression& head, const std::vector<expression>& args) : head_(head), args_(args) {}
		bool is_atom() const { return args_.size() == 0; }

		std::string to_FullForm() const {
			if (is_atom()) {
				return head_.to_FullForm();
			}
			else {
				std::string result = head_.to_FullForm() + "[";
				for (size_t i = 0; i < args_.size(); ++i) {
					result += args_[i].to_FullForm();
					if (i + 1 < args_.size()) {
						result += ", ";
					}
				}
				result += "]";
				return result;
			}
		}
	};

	struct parser {
		lexer lexer_;
		lexer::token currentToken_;

		void consume(lexer::token_type expected) {
			if (currentToken_.type == expected) {
				currentToken_ = lexer_.nextToken();
			}
			else {
				std::cerr << "Parse Error: Unexpected token '" << currentToken_.value
					<< "', expected " << static_cast<int>(expected)
					<< " at position " << currentToken_.position << std::endl;
			}
		}

		atom_expression parse_atom() {
			std::string value = currentToken_.value;

			switch (currentToken_.type) {
			case lexer::IDENTIFIER: {
				consume(lexer::IDENTIFIER);
				return atom_expression(atom_type::Symbol, value);
			}
			case lexer::EXPRESSION: {
				consume(lexer::EXPRESSION);
				return atom_expression(atom_type::Expression, value);
			}
			case lexer::INTEGER: {
				consume(lexer::INTEGER);
				return atom_expression(atom_type::Integer, value);
			}
			case lexer::REAL: {
				consume(lexer::REAL);
				return atom_expression(atom_type::Real, value);
			}
			case lexer::STRING: {
				consume(lexer::STRING);
				return atom_expression(atom_type::String, value);
			}
			default:
				std::cerr << "Parse Error: Unexpected token '" << currentToken_.value
					<< "' at position " << currentToken_.position << std::endl;
				return atom_expression(atom_type::Null, "");
			}
		}

		expression parse_expression() {
			atom_expression head = parse_atom();

			if (currentToken_.type == lexer::LBRACKET) {
				consume(lexer::LBRACKET);
				std::vector<expression> args;

				if (currentToken_.type != lexer::RBRACKET) {
					args.push_back(parse_expression());

					while (currentToken_.type == lexer::COMMA) {
						consume(lexer::COMMA);
						args.push_back(parse_expression());
					}
				}

				// if we meet expr like f[], we use a null atom as argument
				if (args.size() == 0) {
					args.push_back(expression(atom_expression(atom_type::Null, "")));
				}

				consume(lexer::RBRACKET);
				return expression(head, args);
			}

			return head;
		}

		parser(const std::string_view input) : lexer_(input) {
			currentToken_ = lexer_.nextToken();
		}

		expression parse() {
			expression result = parse_expression();

			if (currentToken_.type != lexer::END) {
				std::cerr << "Parse Error: Unexpected token at end: " << currentToken_.value << std::endl;
			}

			return result;
		}
	};

	// we allow use { }, so we need to convert { } to List[ ]
	inline expression parse_FullForm(const std::string_view str) {
		std::string mod_str;
		mod_str.reserve(str.size() + 10);
		for (auto c : str) {
			if (c == '{')
				mod_str += "List[";
			else if (c == '}')
				mod_str += "]";
			else
				mod_str += c;
		}

		parser parser(mod_str);
		return parser.parse();
	}
} // namespace WXF_PARSER::FullForm

namespace WXF_PARSER {
	// we allow use a map to store function that generating sub-expressions
	template<typename Callable>
	void fullform_to_wxf(Encoder& encoder, const FullForm::expression& expr,
		const std::unordered_map<std::string, Callable>& map) {

		if (expr.is_atom()) {
			switch (expr.head_.get_type()) {
			case FullForm::atom_type::Integer:
				encoder.push_integer(std::stoll(expr.head_.get_value()));
				break;
			case FullForm::atom_type::Real:
				encoder.push_real(std::stod(expr.head_.get_value()));
				break;
			case FullForm::atom_type::String:
				encoder.push_string(expr.head_.get_value());
				break;
			case FullForm::atom_type::Symbol:
				encoder.push_symbol(expr.head_.get_value());
				break;
			case FullForm::atom_type::Null:
				break;
			case FullForm::atom_type::Expression: {
				auto& vv = expr.head_.get_value();
				auto it = map.find(vv);
				if (it != map.end())
					it->second(encoder);
				else
					std::cerr << "Error: expression id " << vv << " not found in map." << std::endl;
				break;
			}
			default:
				std::cerr << "Error: unknown atom type." << std::endl;
				break;
			}
		}
		else {
			size_t len = expr.args_.size();
			if (expr.args_[0].is_atom() &&
				expr.args_[0].head_.get_type() == FullForm::atom_type::Null) {
				len = 0;
			}
			std::string_view name = expr.head_.get_value();

			// Rule or DelayRule is not special in FullForm, we treat them as normal functions
			encoder.push_function(name, len);

			for (size_t i = 0; i < len; i++) {
				fullform_to_wxf(encoder, expr.args_[i], map);
			}
		}
	}

	inline void fullform_to_wxf(Encoder& encoder, const FullForm::expression& expr,
		const std::unordered_map<std::string, Encoder>& map) {
		std::unordered_map<std::string, std::function<void(Encoder&)>> func_map;
		for (const auto& [key, value] : map) {
			func_map[key] = [&value](Encoder& enc) {
				enc.push_ustr(value.buffer);
				};
		}
		fullform_to_wxf(encoder, expr, func_map);
	}

	template<typename MapType>
	Encoder fullform_to_wxf(const std::string_view ff_template, const MapType& map, bool include_head = true) {
		Encoder encoder;
		encoder.buffer.reserve(ff_template.size() * 32); // reserve some space
		if (include_head) {
			encoder.buffer.push_back(56); // WXF head
			encoder.buffer.push_back(58); // WXF head
		}
		fullform_to_wxf(encoder, FullForm::parse_FullForm(ff_template), map);
		return encoder;
	}
} // namespace WXF_PARSER