#ifndef __WORD
#define __WORD

#include <unordered_map> // For hash template
#include <ostream> // ostream::<<
#include <iostream> // Debugging
#include <string.h>
#include <math.h>
#include <assert.h>
#include "base_table.h"

template<class Block, int LEN>
class __word
{
	private:
		Block buffer[LEN];
	public:
		
		__word()
		{
			clear();
		};
		
		__word(const std::string &m_seq)
		{
			*this = m_seq;
		};
		
		size_t hash() const; // Specialized in word.cpp
		void push_back(unsigned char m_b); // Specialized in word.cpp
		void pop_front(); // Specialized in word.cpp
		
		// Size of Word intersection
		unsigned int operator&(const __word &m_rhs) const; // Specialized in word.cpp
		
		// Use dynamic programming to find the maximum, ungapped *fractional* overlap between two
		// words
		float max_overlap(const __word &m_rhs) const
		{
			// The "left hand" word is the query and the "right hand" word is 
			// the subject.
			
			// Make sure that scoring matrix can store the highest possible
			// score from two words of the maximum word size.
			typedef unsigned char DPScore;
			
			assert( ( 1 << 8*sizeof(DPScore) ) >= max_size() );
			
			DPScore max_score = 0;
			
			DPScore *dp = new DPScore[max_size()];
			
			if(dp == NULL){
				throw __FILE__ ":max_overlap: Unable to allow dynamic programming row";
			}
			
			// Initialize the dynamic programming row to zero
			memset( dp, 0, max_size() );
			
			const unsigned int query_start = start();
			const unsigned int query_stop = stop();
			
			const unsigned int subject_start = m_rhs.start();
			const unsigned int subject_stop = m_rhs.stop();
			
			for(unsigned int i = query_start;i <= query_stop;++i){
				
				DPScore last_score = 0;
				
				const unsigned char q = get(i);
				
				for(unsigned int j = subject_start;j <= subject_stop;++j){
								
					const DPScore curr_score = dp[j];

					dp[j] = last_score;

					if(m_rhs.get(j) == q){

						++dp[j];

						max_score = (max_score < dp[j]) ? dp[j] : max_score;
					}

					last_score = curr_score;
				}
			}
			
			delete [] dp;
			
			return float(max_score)/std::max( size(), m_rhs.size() );
		};
		
		void shift_left(); // Specialized in word.cpp
		void shift_right(); // Specialized in word.cpp
		
		inline double degeneracy() const
		{
			double ret = 1.0;
			
			for(int i = 0;i < LEN;++i){
				
				const Block &ref = buffer[i];
				
				for(unsigned int j = 0;j < sizeof(Block);++j){
					
					const unsigned char b = (ref >> j*8);
					const unsigned char n1 = b & BASE_MASK;
					const unsigned char n2 = (b >> BITS_PER_BASE) & BASE_MASK;
					
					unsigned char d = 0;
					
					// Count the number of bits in the first nibble (half byte)
					for(int k = 0;k < BITS_PER_BASE;++k){
						d += (n1 >> k) & 1;
					}
					
					// Ignore empty positions
					if(d != 0){
						ret *= d;
					}
					
					d = 0;
					
					// Count the number of bits in the second nibble (half byte)
					for(int k = 0;k < BITS_PER_BASE;++k){
						d += (n2 >> k) & 1;
					}
					
					// Ignore empty positions
					if(d != 0){
						ret *= d;
					}
				}
			}
			
			return ret;
		};
		
		inline __word complement() const
		{
			__word ret;
			
			//      b[0]   b[1]
			// 5'-ABCDEFGabcdefg-3'
			//
			// g = comp(A)
			// A = comp(g)
			// f = comp(B)
			// B = comp(f)
			// ...
			const int first = start();
			const int last = stop();
			
			int dest = 0;
			
			for(int src = last;src >= first;--src, ++dest){
				
				Block comp_b = 0;
				
				const unsigned char b = get(src);
				
				if(b & Base::A){
					comp_b |= Base::T;
				}

				if(b & Base::T){
					comp_b |= Base::A;
				}

				if(b & Base::G){
					comp_b |= Base::C;
				}

				if(b & Base::C){
					comp_b |= Base::G;
				}
				
				ret.set(comp_b, dest);
			}
			
			return ret;
		};
		
		inline bool operator==(const __word &m_rhs) const
		{
			for(int i = 0;i < LEN;++i){
				
				if(buffer[i] != m_rhs.buffer[i]){
					return false;
				}
			}
			
			return true;
		};
		
		inline bool operator<(const __word &m_rhs) const
		{
			for(int i = 0;i < LEN;++i){
				
				if(buffer[i] < m_rhs.buffer[i]){
					return true;
				}
				
				if(buffer[i] > m_rhs.buffer[i]){
					return false;
				}
			}
			
			return false;
		};
		
		inline bool operator>(const __word &m_rhs) const
		{
			for(int i = 0;i < LEN;++i){
				
				if(buffer[i] > m_rhs.buffer[i]){
					return true;
				}
				
				if(buffer[i] < m_rhs.buffer[i]){
					return false;
				}
			}
			
			return false;
		};
		
		inline void clear()
		{
			memset( buffer, 0, LEN*sizeof(Block) );
		};
		
		inline __word& operator=(const std::string &m_seq)
		{
			clear();
			
			const unsigned int len = m_seq.size();
			
			if( len > max_size() ){
				throw __FILE__ ":__word::operator=: Input sequence is too large";
			}
			
			for(unsigned int i = 0;i < len;++i){
				set(base_to_bits(m_seq[i]), i);
			}
			
			return *this;
		};
		
		inline bool empty() const
		{
			return (size() == 0);
		};
		
		inline int start() const
		{
			for(int i = 0;i < LEN;++i){
				
				const Block &src = buffer[i];
				
				for(unsigned int j = 0;j < 8*sizeof(Block);j += BITS_PER_BASE){
					
					if( ( ( src >> (8*sizeof(Block) - BITS_PER_BASE - j) ) & BASE_MASK ) != Base::EOS){
						return i*2*sizeof(Block) + j/BITS_PER_BASE;
					}
				}
			}
			
			return max_size();
		};
		
		inline int stop() const
		{
			for(int i = LEN - 1;i >= 0;--i){
				
				const Block &src = buffer[i];
				
				for(unsigned int j = 0;j < 8*sizeof(Block);j += BITS_PER_BASE){
					
					if( ( (src >> j) & BASE_MASK ) != Base::EOS){
						return i*2*sizeof(Block) + (8*sizeof(Block) - 1 - j)/BITS_PER_BASE;
					}
				}
			}
			
			return -1;
		};
		
		inline unsigned char get(const unsigned int &m_index) const
		{
			assert( m_index < max_size() );
			
			const Block &src = buffer[ m_index/( 2*sizeof(Block) ) ];
			
			return (src >> ( (2*sizeof(Block) - 1) - m_index%( 2*sizeof(Block) ) )*BITS_PER_BASE) & BASE_MASK;
		};
		
		inline std::pair<unsigned char, unsigned char> get_last_two() const
		{
			const int last = stop();
			const int penultimate = last - 1;
			
			return std::make_pair(get(penultimate), get(last) );
		};
		
		inline void set(Block m_base, const unsigned int &m_index)
		{
			assert( m_index < max_size() );
			
			Block &src = buffer[ m_index/( 2*sizeof(Block) ) ];
			const unsigned int offset =  ( (2*sizeof(Block) - 1) - m_index%( 2*sizeof(Block) ) )*BITS_PER_BASE;
			
			src &= ~(Block(Base::N) << offset);
			src |= (m_base << offset);
		};
		
		inline void mask(Block m_base, const unsigned int &m_index)
		{
			assert( m_index < max_size() );
			
			Block &src = buffer[ m_index/( 2*sizeof(Block) ) ];
			const unsigned int offset =  ( (2*sizeof(Block) - 1) - m_index%( 2*sizeof(Block) ) )*BITS_PER_BASE;

			src |= (m_base << offset);
		};
		
		inline void unmask(Block m_base, const unsigned int &m_index)
		{
			assert( m_index < max_size() );
			
			Block &src = buffer[ m_index/( 2*sizeof(Block) ) ];
			const unsigned int offset =  ( (2*sizeof(Block) - 1) - m_index%( 2*sizeof(Block) ) )*BITS_PER_BASE;

			src &= ~(m_base << offset);
		};
		
		inline unsigned char front() const
		{
			return get( start() );
		};
		
		inline unsigned char back() const
		{
			return get( stop() );
		};
		
		unsigned int size() const;
		
		inline unsigned int max_size() const
		{
			return 2*LEN*sizeof(Block);
		};
		
		inline void shrink_front()
		{
			const int i = start();
			
			if( i < int( max_size() ) ){
				set(Base::EOS, i);
			}
		};
		
		inline void shrink_back()
		{
			const int i = stop();
			
			if(i >= 0){
				set(Base::EOS, i);
			}
			
		};
		
		inline void grow_front(unsigned char m_b)
		{
			const int i = start() - 1;
			
			if(i >= 0){
				set(m_b, i);
			}
		};
		
		inline void grow_back(unsigned char m_b)
		{
			const int i = stop() + 1;
			
			if( i < int( max_size() ) ){
				set(m_b, i);
			}
		};
		
		inline void center()
		{
			int left = start();
			int right = stop();
			
			if(left > right){
				// This is an empty word
				return;
			}
			
			right = max_size() - right;
			
			// left == number of gaps on the 5' end
			// right == number of gaps on the 3' end
			const int delta = (right - left)/2;
				
			if(delta > 0){
				for(int i = 0;i < delta;++i){
					shift_right();
				}
			}
			else{
				for(int i = 0;i > delta;--i){
					shift_left();
				}
			}
		};
		
		// Word union
		__word operator|(const __word &m_rhs) const
		{
			__word ret(*this);
			
			const int first = start();
			const int last = stop();
			
			for(int i = first;i <= last;++i){
				
				const unsigned char b = m_rhs.get(i);
				
				if(b == Base::EOS){
					continue;
				}
				
				ret.mask(b, i);
			}
			
			return ret;
		};
		
		// Needs work before it should be used...
		#ifdef INCLUDE_MELTING_TEMPERATURE
		// Return the pair<min tm, max tm> for a perfect match to this oligo
		inline std::pair<float, float> melting_temperature() const
		{
			// For now, use the formula listed below from (http://www.basic.northwestern.edu/biotools/oligocalc.html)
			// Tm = 100.5 + (41.0 * N_GC/N) - (820.0/N) + 16.6*log10([Na+])
			// 1) This formula needs to be verified, because I haven't been able to track down a valid
			// reference for it!
			// 2) Need a simple formula that allows for [strand] correction (due to potentially high degeneracy)
			const float salt_correction = 16.6f*log10f(0.05f);
			
			// The upper bound on the number of GC bases
			unsigned int gc_upper = 0;
			
			// The lower bound on the number of GC bases
			unsigned int gc_lower = 0;
			
			// The total number of ATGC bases
			unsigned int N = 0;
			
			// The oligo degeneracy
			double degen = 1.0;
			
			for(int i = 0;i < LEN;++i){
				
				const Block &ref = buffer[i];
				
				for(unsigned int j = 0;j < sizeof(Block);++j){
					
					const unsigned char b = (ref >> j*8);
					
					// First nibble
					const unsigned char n1 = b & BASE_MASK;
					
					// Second nibble
					const unsigned char n2 = (b >> BITS_PER_BASE) & BASE_MASK;
					
					N += (n1 > 0);
					N += (n2 > 0);
					
					gc_upper += ( (n1 & Base::S) != 0 ) + ( (n2 & Base::S) != 0 ); // S = G | C
					
					// Compute the upper bound on the number of AT bases (which we will invert later
					// to compute the lower bound on the number of GC bases)
					gc_lower += ( (n1 & Base::W) != 0 ) + ( (n2 & Base::W) != 0 ); // W = A | T
					
					unsigned char d = 0;
					
					// Count the number of bits in the first nibble (half byte)
					for(int k = 0;k < BITS_PER_BASE;++k){
						d += (n1 >> k) & 1;
					}
					
					// Ignore empty positions
					if(d != 0){
						degen *= d;
					}
					
					d = 0;
					
					// Count the number of bits in the second nibble (half byte)
					for(int k = 0;k < BITS_PER_BASE;++k){
						d += (n2 >> k) & 1;
					}
					
					// Ignore empty positions
					if(d != 0){
						degen *= d;
					}
				}
			}
			
			gc_lower = N - gc_lower;
						
			return std::make_pair(
					100.5f + (41.0f*gc_lower - 820.0f)/N + salt_correction, 
					100.5f + (41.0f*gc_upper - 820.0f)/N + salt_correction
				);
		};
		#endif // INCLUDE_MELTING_TEMPERATURE
		
		// __word is its own iterator!
		__word begin() const
		{
			__word ret;
			
			for(int i = 0;i < LEN;++i){
				
				const Block &src = buffer[i];
				Block &dst = ret.buffer[i];
				
				for(unsigned int j = 0;j < sizeof(Block);++j){
					
					const unsigned char src_byte = (src >> j*8);
					unsigned char dst_byte = 0;
					
					for(int k = 0;k < BITS_PER_BASE;++k){
						
						const unsigned char curr_base = (1 << k);
						
						if(src_byte & curr_base){
						
							dst_byte |= curr_base;
							break; // Find the first non-zero bit
						}
					}
					
					for(int k = BITS_PER_BASE;k < 8;++k){
						
						const unsigned char curr_base = (1 << k);
						
						if(src_byte & curr_base){
							
							dst_byte |= curr_base;
							break; // Find the first non-zero bit
						}
					}
					
					dst |= ( (Block)dst_byte << j*8);
				}
			}
			
			return ret;
		};
		
		// Return true if the iterator has been set to valid sequence.
		// Return false if the iterator is not valid
		bool next(__word &m_iter) const
		{
			for(int i = 0;i < LEN;++i){
				
				Block &iter = m_iter.buffer[i];
				const Block &src = buffer[i];
				
				for(unsigned int j = 0;j < sizeof(Block);++j){
					
					const unsigned char iter_byte = (iter >> j*8);
					const unsigned char src_byte = (src >> j*8);
					
					// First nibble
					unsigned char iter_nib = iter_byte & BASE_MASK;
					unsigned char src_nib = src_byte & BASE_MASK;
					
					bool wrapped = false;
					
					if(iter_nib != Base::EOS){

						do{
							if(iter_nib == Base::T){

								iter_nib = Base::A;
								wrapped = true;
							}
							else{
								iter_nib <<= 1;
							}

						} while( !(iter_nib & src_nib) );
						
						// Erase the first nibble of the iterator
						iter &= ~(Block(15) << j*8);
						
						// Set the first nibble of the iterator
						iter |= (Block(iter_nib) << j*8);
						
						if(!wrapped){
							return true;
						}
					}
					
					// Second nibble
					iter_nib = (iter_byte >> BITS_PER_BASE) & BASE_MASK;
					src_nib = (src_byte >> BITS_PER_BASE) & BASE_MASK;
					wrapped = false;
					
					if(iter_nib != Base::EOS){

						do{
							if(iter_nib == Base::T){

								iter_nib = Base::A;
								wrapped = true;
							}
							else{
								iter_nib <<= 1;
							}

						} while( !(iter_nib & src_nib) );
						
						// Erase the second nibble of the iterator
						iter &= ~( Block(15) << (j*8 + BITS_PER_BASE) );
						
						// Set the first nibble of the iterator
						iter |= ( Block(iter_nib) << (j*8 + BITS_PER_BASE) );
						
						if(!wrapped){
							return true;
						}
					}
				}
			}
			
			// If we get here, we have return to the begin() iterator and should stop iterating!
			return false;
		};
		
		void write(std::string &m_buffer) const
		{
			const int first = start();
			const int last = stop();

			for(int i = first;i <= last;++i){
				m_buffer.push_back( bits_to_base( get(i) ) );
			}
		}

		inline std::string str() const
		{
			std::string ret;
			
			write(ret);
			
			return ret;
		};
		
		inline size_t mpi_size() const
		{
			return sizeof(Block)*LEN;
		};
		
		inline unsigned char* mpi_pack(unsigned char* m_ptr) const
		{
			memcpy(m_ptr, buffer, sizeof(Block)*LEN);
			m_ptr += sizeof(Block)*LEN;
			
			return m_ptr;
		};
		
		inline unsigned char* mpi_unpack(unsigned char* m_ptr)
		{
			memcpy(buffer, m_ptr, sizeof(Block)*LEN);
			m_ptr += sizeof(Block)*LEN;
			
			return m_ptr;
		};
};

typedef __word<unsigned long int, 2> Word;

#define	WORD_LENGTH	32 // == sizeof(unsigned long int)*2*2 ==  max word capacity

typedef __word<unsigned int, 1> ShortWord;

template<>
size_t __word<unsigned long int, 2>::hash() const;

template<>
size_t __word<unsigned int, 1>::hash() const;

namespace std
{
	template <>
	struct hash < __word<unsigned long int, 2> >
	{
	    size_t operator()(const __word<unsigned long int, 2>& m_word) const
	    {
        	return m_word.hash();
	    }
	};
	
	template <>
	struct hash < __word<unsigned int, 1> >
	{
	    size_t operator()(const __word<unsigned int, 1>& m_word) const
	    {
        	return m_word.hash();
	    }
	};
}

std::ostream& operator << (std::ostream &m_s, const Word &m_w);

void append(std::string &m_buffer, const Word &m_w);

float taq_mama_correction(const std::pair<unsigned char, unsigned char> &m_primer, 
	const std::pair<unsigned char, unsigned char> &m_template_seq);

#endif // __WORD
