
#include <ATEAMS++/ATEAMS++.h>
#include <ATEAMS++/models/Glauber.h>
#include <ATEAMS++/complexes/Cubical.h>
#include <ATEAMS++/statistics/schedules.h>
#include <ATEAMS++/statistics/observables.h>
#include <ATEAMS++/arithmetic/options.h>

using namespace ATEAMS;
using namespace std;

int main() {
	vector<int> corners = {3,3,3,3};
	complexes::Cubical C(corners, true);

	models::GlauberParameters params;
	params.field = 3;
	params.temperatureFunction = statistics::selfdual(params.field);
	params.dimension = 2;
	params.DEBUG = false;

	models::Glauber G(&C, params);
	G.initialize();

	ATEAMS::arithmetic::ThreadOptions options;

	int N = 1000;
	
	for (int i=0; i < N; i++) {
		G.sample(i, options);
		cout << statistics::occupation(G.complex, G.state.cochain, G.field, params.dimension) << endl;
	}
	
	return 0;
}

