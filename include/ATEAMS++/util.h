
#ifndef ATEAMS_UTIL_H
#define ATEAMS_UTIL_H

#include <termios.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <iostream>
#include <thread>

using namespace std;

int getch_key() {
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


void key_listener(std::atomic<bool>& stop_flag) {
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



template <typename T>
void _product(vector<vector<T>> sets, vector<vector<T>> &c, vector<T> &active, int depth) {
	for (int i=0; i<sets[depth].size(); i++) {
		active[depth] = sets[depth][i];
		if (depth < sets.size()-1) _product(sets, c, active, depth+1);
		else c.push_back(active);
	}
}


template <typename T>
vector<vector<T>> product(vector<vector<T>> sets) {
	vector<vector<T>> c;
	vector<T> active(sets.size());
	_product<T>(sets, c, active, 0);

	return c;
}

#endif