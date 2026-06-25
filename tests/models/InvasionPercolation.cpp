
#include "ATEAMS++/ATEAMS++.h"
#include "ATEAMS++/models/InvasionPercolation.h"
#include "ATEAMS++/complexes/Cubical.h"
#include "ATEAMS++/statistics/schedules.h"
#include "ATEAMS++/arithmetic/options.h"

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {3,3,3,3};
	complexes::Cubical C(corners, true);

	models::InvasionPercolationParameters params;
	params.stoppingFunction = arithmetic::stopInvadingAt({3,4});
	params.dimension = 2;

	models::InvasionPercolation P(&C, params);
	ATEAMS::arithmetic::ThreadOptions options;

	int N = 1000;
	vector<double> occupancy(N);
	
	for (int i=0; i < N; i++) {
		P.sample(i, options);
		occupancy[i] = (double)P.state.includes.size()/(double)C.Cells[params.dimension];
	}
	double T = accumulate(occupancy.begin(), occupancy.end(), (double)0);
	double e = (double)T/(double)N;

	cout << format("estimated critical probability is {}", e) << endl;
	
	return 0;
}

