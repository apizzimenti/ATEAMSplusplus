
#include <ATEAMS++/ATEAMS++.h>
#include <ATEAMS++/models/Bernoulli.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <ATEAMS++/statistics/schedules.h>
#include <ATEAMS++/arithmetic/options.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {3,3,3,3};
	complexes::Cubical C(corners, true);

	models::BernoulliParameters params;
	params.p = 0.5;
	params.dimension = 2;
	params.DEBUG = false;

	models::Bernoulli P(&C, params);
	ATEAMS::arithmetic::ThreadOptions options;

	int N = 1000;
	vector<int> rank(N);
	
	for (int i=0; i < N; i++) {
		P.sample(i, options);
		rank[i] = P.state.rank;
	}

	int T = std::accumulate(rank.begin(), rank.end(), 0);
	double e = (double)T/(double)N;

	std::cout << "expected rank is " << std::format("{:.2f}", e) << std::endl;
	
	return 0;
}

