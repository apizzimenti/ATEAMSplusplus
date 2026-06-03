
#ifndef ATEAMS_UTIL_H
#define ATEAMS_UTIL_H

#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <iostream>
#include <thread>

#include "ATEAMS++/common.h"

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
template <typename Matrix>
inline void printSparseRREFmat(Matrix A) {
	int rows, columns;
	rows = A.nrow;
	columns = A.ncol;

	std::vector<std::vector<int>> M(rows, std::vector<int>(columns, 0));

	for (int i=0; i<rows; i++) {
		for (auto [j,v] : A.rows[i]) {
			M[i][j] = (int)v;
		}
	}

	for (int i=0; i<rows; i++) {
		for (int j=0; j<columns; j++) {
			std::cout << M[i][j] << " ";
		}
		std::cout << std::endl;
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
		std::cout << (int)k << ": " << (int)v << ", ";
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