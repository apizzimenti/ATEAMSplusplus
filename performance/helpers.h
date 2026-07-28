
#include <vector>
#include <fstream>
#include <sstream>

// Read in the sample filtrations.
std::vector<std::vector<int>> filtrationData(int size, int SCALE, int DIMENSION, int TRIALS) {
	std::vector<std::vector<int>> filtrations(TRIALS);
	int t = 0;

	std::ifstream filtrationData(std::format("./performance/samples/{}.{}.csv", SCALE, DIMENSION));
	std::string token;
	std::string line;

	while (std::getline(filtrationData,line) && t < TRIALS) {
		std::vector<int> filtration(size);
		std::stringstream ll(line);
		int l = 0;

		while (std::getline(ll, token, ',')) {
			filtration[l] = std::stoi(token);
			l++;
		}

		filtrations[t] = filtration;
		t++;
	}

	return filtrations;
}