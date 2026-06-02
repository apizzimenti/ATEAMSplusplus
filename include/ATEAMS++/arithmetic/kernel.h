
#ifndef ATEAMS_ARITHMETIC_KERNEL_H
#define ATEAMS_ARITHMETIC_KERNEL_H

#include "ATEAMS++/common.h"
#include <vector>
#include <random>
#include <set>

namespace ATEAMS {
	ZpVector randomLinearCombination(
		ZpMatrix kernel,
		const Zp& F,
		std::uniform_int_distribution<int>& intuniform,
		std::mt19937& RNG
	);
	
	ZpVector submatrixKernelSample(
		ZpMatrix coboundary,
		const Zp& F,
		std::set<size_t> exclude,
		std::uniform_int_distribution<int>& intuniform,
		std::mt19937& RNG,
		bool _DEBUG=false
	);
}

#endif
