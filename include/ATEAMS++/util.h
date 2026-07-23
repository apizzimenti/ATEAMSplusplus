
#ifndef ATEAMS_UTIL_H
#define ATEAMS_UTIL_H

#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

#include <iostream>
#include <thread>

#include <fstream>
#include <string>

#include "ATEAMS++/common.h"
#include <SparseRREF/sparse_mat.h>


/** @cond */
inline int getch_key() {
	char buf = 0;
	struct termios old = {};
	if (tcgetattr(STDIN_FILENO, &old) < 0) return -1;

	struct termios new_t = old;
	new_t.c_lflag &= ~(ICANON | ECHO);  
	new_t.c_iflag &= ~IXON;             
	if (tcsetattr(STDIN_FILENO, TCSANOW, &new_t) < 0) return -1;

	int flags = fcntl(STDIN_FILENO, F_GETFL, 0);
	if (flags == -1) return -1;
	fcntl(STDIN_FILENO, F_SETFL, flags | O_NONBLOCK);

	int nread = read(STDIN_FILENO, &buf, 1);

	fcntl(STDIN_FILENO, F_SETFL, flags);
	tcsetattr(STDIN_FILENO, TCSADRAIN, &old);

	if (nread > 0)
		return static_cast<int>(buf);
	else
		return -1; 
}


inline void key_listener(std::atomic<bool>& stop_flag) {
	while (!stop_flag) {
		int c = getch_key();
		if (c == 17) { // Ctrl+Q
			std::cout << "\n[Ctrl+Q] pressed. Stopping...\n";
			stop_flag = true;
			break;
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(50));
	}
}


// Helper method for printing SparseRREF matrices.
template <typename T>
void printSparseRREFmat(ATEAMS::SparseMatrix<T> A) {
	int rows, columns;
	rows = A.nrow;
	columns = A.ncol;

	std::vector<std::vector<int>> M(rows, std::vector<int>(columns, 0));

	std::cout << "   |";
	for (int i=0; i < columns; i++) std::cout << std::format("{:^3}", i);
	std::cout << std::endl;

	for (int i=0; i < columns; i++) std::cout << "———";
	std::cout << std::endl;

	for (int i=0; i<rows; i++) {
		for (auto [j,v] : A.rows[i]) {
			M[i][j] = (int)v;
		}
	}

	for (int i=0; i<rows; i++) {
		std::cout << std::format("{:^3}|", i);
		for (int j=0; j<columns; j++) {
			std::cout << std::format("{:^3}",M[i][j]);
		}
		std::cout << std::endl;
	}
}

// Helper method for printing SparseRREF matrices.
template <typename T>
void printSparseRREFmat(ATEAMS::SparseMatrix<T> A, int brow, bool edited) {
	int rows, columns;
	rows = A.nrow;
	columns = A.ncol;

	std::string col = edited ? "\e[31m" : "\e[34m";

	std::vector<std::vector<int>> M(rows, std::vector<int>(columns, 0));

	std::cout << "   |";
	for (int i=0; i < columns; i++) std::cout << std::format("{:^3}", i);
	std::cout << std::endl;

	for (int i=0; i < columns+1; i++) std::cout << "———";
	std::cout << std::endl;

	for (int i=0; i<rows; i++) {
		for (auto [j,v] : A.rows[i]) {
			M[i][j] = (int)v;
		}
	}

	for (int i=0; i<rows; i++) {
		if (i == brow) std::cout << col << "\e[1m";
		std::cout << std::format("{:^3}|", i);
		for (int j=0; j<columns; j++) {
			std::cout << std::format("{:^3}",M[i][j]);
		}
		std::cout << "\e[0m\e[0m" << std::endl;
	}
}

// Helper method for printing SparseRREF matrices.
template <typename T>
void printSparseRREFmat(ATEAMS::SparseMatrix<T> A, int fixed, int edited) {
	int rows, columns;
	rows = A.nrow;
	columns = A.ncol;

	std::string fixedcol = "\e[31m";
	std::string editedcol = "\e[32m";

	std::vector<std::vector<int>> M(rows, std::vector<int>(columns, 0));

	std::cout << "   |";
	for (int i=0; i < columns; i++) std::cout << std::format("{:^3}", i);
	std::cout << std::endl;

	for (int i=0; i < columns+1; i++) std::cout << "———";
	std::cout << std::endl;

	for (int i=0; i<rows; i++) {
		for (auto [j,v] : A.rows[i]) {
			M[i][j] = (int)v;
		}
	}

	for (int i=0; i<rows; i++) {
		if (i == fixed) std::cout << fixedcol << "\e[1m";
		if (i == edited) std::cout << editedcol << "\e[1m";

		std::cout << std::format("{:^3}|", i);
		for (int j=0; j<columns; j++) {
			std::cout << std::format("{:^3}",M[i][j]);
		}

		if (i == fixed) std::cout << "<<< youngest \e[0m";
		if (i == edited) std::cout << "<<< cell \e[0m";
		std::cout << "\e[0m" << std::endl;
	}
}


// Helper method for printing SparseRREF matrices.
template <typename T>
void printSparseRREFmat(ATEAMS::SparseMatrix<T> A, std::vector<int> nextColumnAdded, int zeroed) {
	int rows, columns;
	rows = A.nrow;
	columns = A.ncol;

	std::string zeroedcol = "\e[34m";

	std::vector<std::vector<int>> M(rows, std::vector<int>(columns, 0));

	std::cout << "   |";
	for (int i=0; i < columns; i++) std::cout << std::format("{:^3}", i);
	std::cout << std::endl;

	for (int i=0; i < columns+1; i++) std::cout << "———";
	std::cout << std::endl;

	for (int i=0; i<rows; i++) {
		for (auto [j,v] : A.rows[i]) {
			M[i][j] = (int)v;
		}
	}

	for (int i=0; i<rows; i++) {
		if (i == zeroed) std::cout << zeroedcol << "\e[1m";
		if (i == nextColumnAdded[zeroed]) std::cout << "\e[31m";

		std::cout << std::format("{:^3}|", i);
		for (int j=0; j<columns; j++) {
			std::cout << std::format("{:^3}",M[i][j]);
		}

		std::cout << std::format(" | {:^3}", nextColumnAdded[i]);

		if (i == nextColumnAdded[zeroed]) std::cout << "<<< next added \e[0m";
		if (i == zeroed) std::cout << "<<< zeroed \e[0m";
		std::cout << "\e[0m" << std::endl;
	}
}


// Utility functions for suppressing output to stderr (specifically because SpaSM
// outputs diagnostic information at every step!).
inline int _suppress() {
	fflush(stderr);
	int fd = dup(STDERR_FILENO);
	freopen("/dev/null", "w", stderr);
	return fd;
}

inline void _resume(int fd) {
	fflush(stderr);
	dup2(fd, fileno(stderr));
	close(fd);
}


// Helper method for printing `map`s.'
template <typename MapStorage>
inline void printmap(MapStorage m) {
	std::cout << "{ ";
	for (const auto& [k, v] : m) {
		std::cout << k << ": " << v << ", ";
	}
	std::cout << " }" << std::endl;
}

// Helper method for printing `set`s.'
template <typename t>
void printset(std::set<t> A) {
	std::cout << "{ ";
	for (const auto& a : A) {
		std::cout << (int)a << ", ";
	}
	std::cout << " }" << std::endl;
}

// Helper method for printing `vector`s.
template <typename t>
inline void printvector(std::vector<t> v) {
	std::cout << "[ ";
	for (auto &k : v) {
		std::cout << k << " ";
	}
	std::cout << "]" << std::endl;
}


template <typename MapStorage>
inline void printvectormap(MapStorage m) {
	std::cout << "{ ";
	for (const auto& [k, v] : m) {
		std::cout << "[ ";
		for (auto &p : k) {
			std::cout << p << " ";
		}
		std::cout << "]";
		
		std::cout << ": " << (int)v << ", ";
	}
	std::cout << " }" << std::endl;
}

/** @endcond */


#endif