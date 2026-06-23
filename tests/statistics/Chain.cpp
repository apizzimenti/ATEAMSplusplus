
#include <ATEAMS++/models/Bernoulli.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <ATEAMS++/statistics/Chain.h>

#include <chrono>

using namespace ATEAMS;
using namespace std;

int main() {
	// Construct a cubical complex.
	vector<int> corners = {6,6,6,6};
	complexes::Cubical C(corners, true);

	// Parametrize the model.
	models::BernoulliParameters params;
	params.dimension = 2;

	// Construct the model, then the chain.
	models::Bernoulli percolation(&C, params);

	int iterations = 1000;

	using Chain = statistics::Chain<models::Bernoulli>;
	Chain M(&percolation, iterations);

	// Iterate over the chain, storing the rank of the 2nd homology group as we
	// go.
	vector<int> ranks;

	auto start = chrono::high_resolution_clock::now();

	for (models::BernoulliState* state : M.simulate<models::BernoulliState>()) {
		ranks.push_back(state->rank);
	}

	auto end = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(end-start);
	std::cout << (double)duration.count()/(double)1000 << endl;

	// Compute the expected rank of the 2nd homology group.
	int T = std::accumulate(ranks.begin(), ranks.end(), 0);
	double e = (double)T/(double)iterations;

	std::cout << "expected rank is " << std::format("{:.2f}", e) << std::endl;
	
	return 0;
}


