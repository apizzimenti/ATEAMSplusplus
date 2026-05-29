
#include "ATEAMS++/complexes/Cubical.h"
#include <iostream>

using namespace std;

int main() {
	vector<int> corners = {3,3,3};
	Cubical C(corners, true);
	cout << C.kind << endl;
	cout << C.periodic << endl;
	return 0;
}
