
#include <ATEAMS++/ATEAMS++.h>
#include <ATEAMS++/models/Bernoulli.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <ATEAMS++/statistics/schedules.h>
#include <ATEAMS++/arithmetic/options.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {10,10,10,10};
	complexes::Cubical C(corners, true);

	models::BernoulliParameters params;
	params.p = 0.5;
	params.dimension = 1;
	params.DEBUG = true;

	models::Bernoulli P(&C, params);
	ATEAMS::arithmetic::ThreadOptions options;
	
	for (int i=0; i < 1000; i++) {
		vector<int> r = P.sample(i, options);
	}
	
	return 0;
}

