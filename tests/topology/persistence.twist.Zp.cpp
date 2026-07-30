
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;

template <typename RingLike>
vector<int> restrictedWrapper(
	ATEAMS::complexes::Complex<RingLike>* complex,
	std::vector<int>& filtration,
	ATEAMS::Ring* R,
	int dimension,
	ATEAMS::arithmetic::ComputeResources<RingLike>& options
) {
	return topology::persistence::twist<RingLike>(
		complex,
		filtration,
		R,
		dimension,
		options
	);
}

template <typename RingLike>
vector<int> fullWrapper(
	ATEAMS::complexes::Complex<RingLike>* complex,
	std::vector<int>& filtration,
	ATEAMS::Ring* R,
	ATEAMS::arithmetic::ComputeResources<RingLike>& options
) {
	return topology::persistence::twist<RingLike>(
		complex,
		filtration,
		R,
		options
	);
}


int main(int argc, char *argv[]) {
	using R = Zp;
	return persistenceDispatcher<R>(argc, argv, restrictedWrapper<R>, fullWrapper<R>);
}