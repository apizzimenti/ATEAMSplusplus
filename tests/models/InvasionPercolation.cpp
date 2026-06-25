
#include "ATEAMS++/ATEAMS++.h"
#include "ATEAMS++/models/InvasionPercolation.h"
#include "ATEAMS++/complexes/Cubical.h"
#include "ATEAMS++/statistics/schedules.h"
#include "ATEAMS++/arithmetic/options.h"

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {3,3};
	complexes::Cubical C(corners, true);

	models::InvasionPercolationParameters params;
	params.stoppingFunction = arithmetic::stopInvadingAt({1});
	params.dimension = 1;

	models::InvasionPercolation P(&C, params);
	ATEAMS::arithmetic::ThreadOptions options;

	int N = 2;
	
	for (int i=0; i < N; i++) {
		P.sample(i, options);
	}
	
	return 0;
}

