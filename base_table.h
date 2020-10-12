#ifndef __BASE_TABLE
#define  __BASE_TABLE

#include <iostream>

#define	BITS_PER_BASE	4
#define	BASE_MASK	0xF

namespace Base
{
	enum {
		EOS = 0, // End-of-sequence
		A = (1 << 0), 
		C = (1 << 1), 
		G = (1 << 2), 
		T = (1 << 3), 
		M = (A | C),
		R = (G | A),
		S = (G | C),
		V = (G | C | A),
		W = (A | T),
		Y = (T | C),
		H = (A | C | T),
		K = (G | T),
		D = (G | A | T),
		B = (G | T | C),
		N = (A | T | C | G)
	};
}

inline unsigned char base_to_bits(char m_base)
{
	switch(m_base){
		case 'A': case 'a':
			return Base::A;
		case 'T': case 't':
		case 'U': case 'u':
			return Base::T;
		case 'G': case 'g':
			return Base::G;
		case 'C': case 'c':
			return Base::C;
		case 'M': case 'm':
			return Base::M;
		case 'R': case 'r':
			return Base::R;
		case 'S': case 's':
			return Base::S;
		case 'V': case 'v':
			return Base::V;
		case 'W': case 'w':
			return Base::W;
		case 'Y': case 'y':
			return Base::Y;
		case 'H': case 'h':
			return Base::H;
		case 'K': case 'k':
			return Base::K;
		case 'D': case 'd':
			return Base::D;
		case 'B': case 'b':
			return Base::B;
		case 'N': case 'n':
		case 'I': case 'i': // For now, treat inosine as an 'N'
		case 'X': case 'x':
			return Base::N;
		case '-':
			return Base::EOS;
		default:
			std::cerr << "base_to_bits() does not understand the symbol \'" 
				<< m_base << "\' (" << int(m_base) << ")" << std::endl;

			throw __FILE__ ":base_to_bits: Illegal base";
			break;
	};
};

inline char bits_to_base(unsigned char m_base)
{
	switch(m_base){
		case Base::A:
			return 'A';
		case Base::T:
			return 'T';
		case Base::G:
			return 'G';
		case Base::C:
			return 'C';
		case Base::M:
			return 'M';
		case Base::R:
			return 'R';
		case Base::S:
			return 'S';
		case Base::V:
			return 'V';
		case Base::W:
			return 'W';
		case Base::Y:
			return 'Y';
		case Base::H:
			return 'H';
		case Base::K:
			return 'K';
		case Base::D:
			return 'D';
		case Base::B:
			return 'B';
		case Base::N:
			return 'N';
		case Base::EOS:
			return '-';
		default:
			std::cerr << "m_base = " << int(m_base) << std::endl;
			
			throw __FILE__ ":bits_to_base: Illegal base";
			break;
	};

	// We should never get here
	return '?';
};

inline bool is_degen(unsigned char m_base)
{
	switch(m_base){
		case Base::A:
		case Base::T:
		case Base::G:
		case Base::C:
			return false;
		default:
			return true;
	};
	
	return true;
}

#endif // __BASE_TABLE
