
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
	ATEAMS::arithmetic::ComputeResources<RingLike>& resources
) {
	return topology::persistence::PHAT<RingLike>(
		complex,
		filtration,
		dimension,
		resources
	);
}

template <typename RingLike>
vector<int> fullWrapper(
	ATEAMS::complexes::Complex<RingLike>* complex,
	std::vector<int>& filtration,
	ATEAMS::Ring* R,
	ATEAMS::arithmetic::ComputeResources<RingLike>& resources
) {
	return topology::persistence::PHAT<RingLike>(
		complex,
		filtration,
		resources
	);
}


int main(int argc, char *argv[]) {
	using R = Zp;
	return persistenceDispatcher<R>(argc, argv, restrictedWrapper<R>, fullWrapper<R>);
}