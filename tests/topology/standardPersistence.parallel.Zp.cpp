
#include <ATEAMS++/ATEAMS++.h>
#include "tests.h"

using namespace ATEAMS;
using namespace std;


int main(int argc, char *argv[]) {
	using R = Zp;
	return persistenceDispatcher<R>(argc, argv, topology::standardPersistence<R>);
}