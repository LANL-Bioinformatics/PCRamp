#include "mpi_util.h"
#include "pcramp.h"
#include "bitset.h"
#include "sequence.h"
#include "assay.h"
#include <string.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for std::string
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<string>(const string &m_str)
{
	return sizeof(size_t) + m_str.size();
}

template<>
unsigned char* mpi_pack<string>(unsigned char* m_ptr, const string &m_str)
{
	size_t len = m_str.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	memcpy(m_ptr, m_str.c_str(), len);
	m_ptr += len;
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<string>(unsigned char* m_ptr, string &m_str)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_str.assign( (char*)m_ptr, len );
	m_ptr += len;
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for __uint128_t
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<__uint128_t>(const __uint128_t &m_obj)
{
	return sizeof(__uint128_t);
}

template<>
unsigned char* mpi_pack<__uint128_t>(unsigned char* m_ptr, const __uint128_t &m_obj)
{
	memcpy( m_ptr, &m_obj, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<__uint128_t>(unsigned char* m_ptr, __uint128_t &m_obj)
{
	memcpy( &m_obj, m_ptr, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Options
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const Options &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		OPTIONS_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, Options &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		OPTIONS_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const Options &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		OPTIONS_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Score
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const Score &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		SCORE_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, Score &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		SCORE_MEMBERS
        #undef VARIABLE
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const Score &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		SCORE_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for BitSet
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const BitSet &m_obj)
{
	const size_t num_bits = m_obj.size();
	
	return sizeof(size_t) + num_bits/8 + ( (num_bits%8 > 0) ? 1: 0 );
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, BitSet &m_obj)
{
	size_t num_bits = 0;
			
	memcpy(&num_bits,  m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);

	m_obj.resize(num_bits);

	unsigned char current_bit = 7; // Start at the left-most bit
	
	for(size_t i = 0;i < num_bits;++i){

		m_obj[i] = (*m_ptr >> current_bit) & 1;
		
		if(current_bit == 0){
		
			++m_ptr;
			current_bit = 7;
		}
		else{
			--current_bit;
		}
	}

	// mpi_pack and mpi_unpack are byte aligned, so skip any
	// remaining bits
	if(current_bit != 7){
		++m_ptr;
	}
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const BitSet &m_obj)
{

	size_t num_bits = m_obj.size();
			
	memcpy( m_ptr, &num_bits, sizeof(size_t) );
	m_ptr += sizeof(size_t);

	unsigned char current_bit = 7; // Start at the left-most bit
	
	for(size_t i = 0;i < num_bits;++i){

		// Before setting the left-most bit, we need
		// to initialize the entire byte to zero
		if(current_bit == 7){
			*m_ptr = 0;
		}
	
		*m_ptr |= (m_obj[i] << current_bit);
		
		if(current_bit == 0){
		
			++m_ptr;
			current_bit = 7;
		}
		else{
			--current_bit;
		}
	}

	// mpi_pack and mpi_unpack are byte aligned, so skip any
	// remaining bits
	if(current_bit != 7){
		++m_ptr;
	}
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Sequence
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const Sequence &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		SEQUENCE_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, Sequence &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		SEQUENCE_MEMBERS
        #undef VARIABLE

	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const Sequence &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		SEQUENCE_MEMBERS
	#undef VARIABLE

	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for AmpliconBounds
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const AmpliconBounds &m_obj)
{
	size_t ret = 0;

	#define VARIABLE(A, B) ret += mpi_size(m_obj.B);
		AMPLICON_BOUNDS_MEMBERS
	#undef VARIABLE

	return ret;	
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, AmpliconBounds &m_obj)
{
	#define VARIABLE(A, B) m_ptr = mpi_unpack(m_ptr, m_obj.B);
		AMPLICON_BOUNDS_MEMBERS
    #undef VARIABLE

	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const AmpliconBounds &m_obj)
{

	#define VARIABLE(A, B) m_ptr = mpi_pack(m_ptr, m_obj.B);
		AMPLICON_BOUNDS_MEMBERS
	#undef VARIABLE

	return m_ptr;
}
