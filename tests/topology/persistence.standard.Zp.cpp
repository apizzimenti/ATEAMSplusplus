
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

template <typename RingLike>
vector<int> wrapper(
	ATEAMS::complexes::Complex<RingLike>* complex,
	std::vector<int>& filtration,
	ATEAMS::Ring* R,
	int dimension,
	ATEAMS::arithmetic::ComputeOptions<RingLike>& options
) {
	return topology::persistence::standard<RingLike>(
		complex,
		filtration,
		R,
		dimension,
		options
	);
}


int main(int argc, char *argv[]) {
	using R = Zp;
	return persistenceDispatcher<R>(argc, argv, wrapper<R>);
}