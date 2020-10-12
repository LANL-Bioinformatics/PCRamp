#include "word.h"
#include <emmintrin.h>

#include <iostream>

using namespace std;

#ifdef POPCNT
	#define	POPCOUNT(X)	__builtin_popcountl(X)
#else
	#define	POPCOUNT(X)	local_popcnt(X)
#endif // POPCNT

// Data layout:
//
//      buffer[0]          buffer[1] ...
// MSB <--------->LSB MSB <--------->LSB ...
// 5'------------------------------------...3'
// front---------------------------------...back

/////////////////////////////////////////////////////////////////////////////
// Specializations for Block == unsigned long int
/////////////////////////////////////////////////////////////////////////////
template<>
size_t __word<unsigned long int, 2>::hash() const
{
	return buffer[0] ^ buffer[1];
}


template<>
void __word<unsigned long int, 2>::push_back(unsigned char m_b)
{
	// Are we still filling the buffer?
	const int last = stop() + 1;
	
	if( last < int( max_size() ) ){
	
		set(m_b, last);
		return;
	}
	
	const unsigned char carry = (buffer[1] >> 60) & BASE_MASK;
	
	// Shift to the left and add at the back
	buffer[1] = (buffer[1] << BITS_PER_BASE) | m_b;
	buffer[0] = (buffer[0] << BITS_PER_BASE) | carry;	
}

// From Wikipedia (popcount_3): http://en.wikipedia.org/wiki/Hamming_weight
inline unsigned int local_popcnt(unsigned long int m_buffer)
{
	m_buffer -= (m_buffer >> 1) & 0x5555555555555555;
	m_buffer = (m_buffer & 0x3333333333333333) + ((m_buffer >> 2) & 0x3333333333333333);
	m_buffer = (m_buffer + (m_buffer >> 4)) & 0x0f0f0f0f0f0f0f0f;
	return (m_buffer * 0x0101010101010101) >> 56;
}

#ifdef SSE_VERSION
union sse_elem {
	
	__m128i sse;
	unsigned long int raw[2];
};
#endif // SSE_VERSION

template<>
unsigned int __word<unsigned long int, 2>::operator&(const __word<unsigned long int, 2> &m_rhs) const
{	 

	#ifdef SSE_VERSION
	// This hand-coded sse version is not any faster than the non-sse version shown
	// below.
	const __m128i zero = _mm_setzero_si128();
	const __m128i mask_low =  _mm_set1_epi16(0x0F0F);
	const __m128i mask_high = _mm_set1_epi16(0xF0F0);

	__m128i a;
	sse_elem b;
	
	// May need to use _mm_loadu_si128 for unaligned load
	a = _mm_loadu_si128( (__m128i*)buffer );
	b.sse = _mm_loadu_si128( (__m128i*)m_rhs.buffer );
	
	// Overlap between words
	a = _mm_and_si128(a, b.sse);
	
	// Problem -- we can't use _mm_cmpgt_epi8 to compare nibble values to zero, since
	// the sse comparison instructions work on *signed* 8-bit values
	// Solution -- compare equal to zero and subtract the number of bits in the result
	// from result we would expect if *all* bits were zero.
	
	// Test the first (low-order) nibble of each byte
	b.sse = _mm_cmpeq_epi8( _mm_and_si128(a, mask_low), zero);
	
	unsigned int count = POPCOUNT(b.raw[0]) + __builtin_popcountl(b.raw[1]);
	
	// Test the second (high-order) nibble of each byte
	b.sse = _mm_cmpeq_epi8( _mm_and_si128(a, mask_high), zero);
	
	count += POPCOUNT(b.raw[0]) + __builtin_popcountl(b.raw[1]);
	
	// Count contains the number of bits that are zero (we are interested in the
	// number of bits that *aren't* zero). Divide the count by 8 to convert from bit
	// count to base count.
	count = (256 - count) >> 3;
	
	return count;
	#endif // SSE_VERSION
	
	const unsigned long int A = buffer[0] & m_rhs.buffer[0];
	const unsigned long int B = buffer[1] & m_rhs.buffer[1];
	
	#ifdef VERSION1
	unsigned long int C =
		     (A & 0x8888888888888888);      // 1000:1000:1000:...
		C |= (A & 0x4444444444444444) << 1; // 0100:0100:0100:...
		C |= (A & 0x2222222222222222) << 2; // 0010:0010:0010:...
		C |= (A & 0x1111111111111111) << 3; // 0001:0001:0001:...

		C |= (B & 0x8888888888888888) >> 1; // 1000:1000:1000:...
		C |= (B & 0x4444444444444444);      // 0100:0100:0100:...
		C |= (B & 0x2222222222222222) << 1; // 0010:0010:0010:...
		C |= (B & 0x1111111111111111) << 2; // 0001:0001:0001:...
	
	return POPCOUNT(C);
	#endif // VERSION1
	
	// This is a "builtin" (GCC only) function that wraps the
	// assembly instruction popcnt or falls back to a library call for
	// counting the number of 1's in a binary number. The assembly 
	// instruction is *much* faster than my original implementation, but is
	// only supported on Nehalem (Intel)/Barcelona (AMD) chips (and newer).
	// Version 2 (shown below) is about the same speed as version 1 (but more
	// compact/harder to read).
	
	// #ifdef POPCNT
// 	return __builtin_popcountl(
// 	#else
// 	return local_popcnt(
// 	#endif // POPCNT
// 		    ((B & 0x8888888888888888) >> 1) |
// 		     (A & 0x8888888888888888) |
// 		     (B & 0x4444444444444444) |
// 		   (((A & 0x4444444444444444) | (B & 0x2222222222222222)) << 1) |
// 		   (((A & 0x2222222222222222) | (B & 0x1111111111111111)) << 2) | 
// 		    ((A & 0x1111111111111111) << 3)
// 	);
	
	// A little faster than the version above (only applies two bitmasks instead of eight).
	return POPCOUNT(
		    ( ( A | (A << 1) | (A << 2) | (A << 3) ) & 0x8888888888888888 ) |  // 1000:1000:1000:...
		    ( ( (B >> 1) | B | (B << 1) | (B << 2) ) & 0x4444444444444444 )    // 0100:0100:0100:...
	);
	
	#ifdef ORIGINAL
	const unsigned long int A = buffer[0] & m_rhs.buffer[0];
	const unsigned long int B = buffer[1] & m_rhs.buffer[1];

	return 
		( (A & 0x000000000000000F) != 0 ) + 
		( (A & 0x00000000000000F0) != 0 ) + 
		( (A & 0x0000000000000F00) != 0 ) + 
		( (A & 0x000000000000F000) != 0 ) + 
		( (A & 0x00000000000F0000) != 0 ) + 
		( (A & 0x0000000000F00000) != 0 ) + 
		( (A & 0x000000000F000000) != 0 ) + 
		( (A & 0x00000000F0000000) != 0 ) + 
		( (A & 0x0000000F00000000) != 0 ) + 
		( (A & 0x000000F000000000) != 0 ) + 
		( (A & 0x00000F0000000000) != 0 ) + 
		( (A & 0x0000F00000000000) != 0 ) + 
		( (A & 0x000F000000000000) != 0 ) + 
		( (A & 0x00F0000000000000) != 0 ) + 
		( (A & 0x0F00000000000000) != 0 ) + 
		( (A & 0xF000000000000000) != 0 ) +
		
		( (B & 0x000000000000000F) != 0 ) + 
		( (B & 0x00000000000000F0) != 0 ) + 
		( (B & 0x0000000000000F00) != 0 ) + 
		( (B & 0x000000000000F000) != 0 ) + 
		( (B & 0x00000000000F0000) != 0 ) + 
		( (B & 0x0000000000F00000) != 0 ) + 
		( (B & 0x000000000F000000) != 0 ) + 
		( (B & 0x00000000F0000000) != 0 ) + 
		( (B & 0x0000000F00000000) != 0 ) + 
		( (B & 0x000000F000000000) != 0 ) + 
		( (B & 0x00000F0000000000) != 0 ) + 
		( (B & 0x0000F00000000000) != 0 ) + 
		( (B & 0x000F000000000000) != 0 ) + 
		( (B & 0x00F0000000000000) != 0 ) + 
		( (B & 0x0F00000000000000) != 0 ) + 
		( (B & 0xF000000000000000) != 0 );
	
	#endif// ORIGINAL
}

template<>
unsigned int __word<unsigned long int, 2>::size() const
{
	
	// Compute the size of the word by masking bits to identify the base positions with
	// one or more bits that have been set. Bases that are equal to Base::EOS (i.e. zero) will
	// be the only bases that have *no* bits set. The masking works by collapsing the four bits
	// corresponding to a single base into a single bit (using successive binary OR operations).
	const unsigned long int A = 
		  ( ( buffer[0] | (buffer[0] << 2) ) & 0xCCCCCCCCCCCCCCCC )  // 1100:1100:1100:...
		| ( ( buffer[1] | (buffer[1] >> 2) ) & 0x3333333333333333 ); // 0011:0011:0011:...
	
	return POPCOUNT(
			(A | ( A << 1) ) & 0xAAAAAAAAAAAAAAAA //1010:1010:1010:...
	);
}

template<>
void __word<unsigned long int, 2>::shift_left()
{
	const unsigned char carry = (buffer[1] >> 60) & BASE_MASK;
	
	buffer[0] = (buffer[0] << BITS_PER_BASE) | carry;
	buffer[1] <<= BITS_PER_BASE;
}

template<>
void __word<unsigned long int, 2>::shift_right()
{
	const unsigned long int carry = buffer[0] & BASE_MASK;
	
	buffer[0] >>= BITS_PER_BASE;
	buffer[1] = (buffer[1] >> BITS_PER_BASE) | (carry << 60);
}

inline int taq_mama_index(unsigned char m_b)
{
	switch(m_b){
		case Base::C:
			return 0;
		case Base::G:
			return 1;
		case Base::A:
			return 2;
		case Base::T:
			return 3;
	};
	
	return -1;
}

float taq_mama_correction(const std::pair<unsigned char, unsigned char> &m_primer, 
	const std::pair<unsigned char, unsigned char> &m_template_seq)
{
	const unsigned int N = 16;
	
	// Construct a look-up table copied from 
	// Table 2 in Li et al. Genomics 83 (2004) 311-320
	// 1) Primer sequence is in columns, template sequence is in rows
	// 2) The Li et. al. paper using the following order for both columns and
	// rows: {CC, GC, AC, TC, CG, GG, AG, TG, CA, GA, AA, TA, CT, GT, AT, TT}
	const float TaqMAMA[] = {
		1.000f, 0.968f, 0.947f, 1.034f, 0.547f, 0.253f, 0.230f, 0.359f, 0.606f, 0.282f, 0.372f, 0.347f, 0.957f, 0.382f, 0.399f, 0.687f, 
		0.989f, 1.000f, 1.023f, 1.000f, 0.420f, 0.662f, 0.445f, 0.367f, 0.870f, 0.512f, 0.492f, 0.508f, 0.372f, 1.000f, 0.492f, 0.714f, 
		1.011f, 1.000f, 1.000f, 1.000f, 0.459f, 0.277f, 0.570f, 0.343f, 0.927f, 0.362f, 0.590f, 0.542f, 0.439f, 0.488f, 0.978f, 0.662f, 
		1.000f, 0.907f, 1.000f, 1.000f, 0.382f, 0.234f, 0.228f, 0.542f, 0.763f, 0.309f, 0.410f, 0.473f, 0.426f, 0.347f, 0.423f, 0.947f, 
		0.590f, 0.334f, 0.445f, 0.323f, 1.000f, 0.978f, 0.927f, 0.989f, 0.907f, 0.645f, 0.525f, 0.455f, 0.927f, 0.408f, 0.408f, 0.707f, 
		0.327f, 0.595f, 0.319f, 0.396f, 0.947f, 1.000f, 0.978f, 0.989f, 0.405f, 0.861f, 0.681f, 0.512f, 0.410f, 0.968f, 0.452f, 0.714f, 
		0.410f, 0.420f, 0.590f, 0.311f, 1.023f, 1.000f, 1.000f, 1.000f, 0.488f, 0.898f, 0.907f, 0.566f, 0.442f, 0.449f, 0.989f, 0.707f, 
		0.423f, 0.343f, 0.305f, 0.585f, 1.034f, 0.879f, 0.927f, 1.000f, 0.473f, 0.720f, 0.547f, 0.957f, 0.459f, 0.374f, 0.459f, 1.023f, 
		1.023f, 0.429f, 0.473f, 0.477f, 1.023f, 0.466f, 0.420f, 0.477f, 1.000f, 0.978f, 0.907f, 0.978f, 0.907f, 0.380f, 0.525f, 0.669f, 
		0.442f, 1.046f, 0.455f, 0.470f, 0.432f, 1.058f, 0.481f, 0.485f, 0.917f, 1.000f, 1.023f, 1.023f, 0.336f, 0.968f, 0.534f, 0.639f, 
		0.617f, 0.452f, 1.011f, 0.439f, 0.492f, 0.504f, 0.978f, 0.462f, 0.989f, 0.947f, 1.000f, 0.978f, 0.405f, 0.405f, 0.888f, 0.606f, 
		0.601f, 0.377f, 0.377f, 1.046f, 0.500f, 0.399f, 0.408f, 1.034f, 0.978f, 0.720f, 0.870f, 1.000f, 0.402f, 0.313f, 0.651f, 0.927f, 
		0.978f, 0.462f, 0.466f, 0.488f, 0.420f, 0.239f, 0.225f, 0.336f, 0.504f, 0.269f, 0.319f, 0.656f, 1.000f, 0.835f, 0.907f, 1.034f, 
		0.429f, 1.011f, 0.473f, 0.477f, 0.340f, 0.413f, 0.357f, 0.354f, 0.352f, 0.538f, 0.413f, 0.794f, 0.927f, 1.000f, 1.058f, 1.000f, 
		0.595f, 0.492f, 0.968f, 0.485f, 0.367f, 0.282f, 0.388f, 0.439f, 0.413f, 0.309f, 0.566f, 0.917f, 0.957f, 0.957f, 1.000f, 0.989f, 
		0.590f, 0.380f, 0.410f, 0.968f, 0.364f, 0.223f, 0.230f, 0.416f, 0.321f, 0.239f, 0.301f, 0.645f, 0.978f, 0.714f, 0.947f, 1.000f

	};
	
	const int a = taq_mama_index(m_primer.first);
	const int b = taq_mama_index(m_primer.second);
	const int c = taq_mama_index(m_template_seq.first);
	const int d = taq_mama_index(m_template_seq.second);
	
	if( (a < 0) || (b < 0) || (c < 0) || (d < 0) ){
		// One or more of the input bases are degenerate -- no correction is applied
		return 1.0f;
	}
	
	const int index_primer = 4*b + a;
	const int index_template = 4*d + c;

	// Clamp the correction at 1.0
	return std::min(1.0f, TaqMAMA[N*index_template + index_primer]);
}

std::ostream& operator << (std::ostream &m_s, const Word &m_w)
{
	const int first =  m_w.start();
	const int last = m_w.stop();
	
	//#define PAD_WORDS
	
	#ifdef PAD_WORDS
	for(int i = 0;i < first;++i){
		m_s << '-';
	}
	#endif // PAD_WORDS
	
	for(int i = first;i <= last;++i){
		m_s << bits_to_base( m_w.get(i) );
	}
	
	#ifdef PAD_WORDS
	for(int i = last + 1;i < int( m_w.max_size() );++i){
		m_s << '-';
	}
	#endif // PAD_WORDS
	
	return m_s;
}
