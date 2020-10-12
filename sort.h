#ifndef __SORT
#define __SORT

// Note that on OS X, the clang compiler defines __GNUC__
#if defined(_OPENMP) && !defined(__clang__)
	#include <parallel/algorithm>

	// Enable OpenMP-based parallel sorting	
	#define	SORT	__gnu_parallel::sort
#else
	#include <algorithm>
	
	// Use standard serial-based sorting
	#define	SORT	std::sort
#endif // _OPENMP

#endif // __SORT
