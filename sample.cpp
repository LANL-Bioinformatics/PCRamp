#include "pcramp.h"

using namespace std;

// Generate a random sequence location in the interval [m_start, m_stop]
int random_location(const int &m_start, const int &m_stop, unsigned int *m_seed_ptr)
{
	assert( (m_start <= m_stop) && m_seed_ptr );

	// Unbiased sampling
	return m_start + rand_r(m_seed_ptr)%(m_stop - m_start);
}

