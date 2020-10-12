#ifndef __SEQUENCE
#define __SEQUENCE

#include <vector>
#include <unordered_map>
#include <string>
#include "base_table.h"
#include "word.h"
#include "mpi_util.h"

#include "read_only_multimap.h"

#define MAP		std::unordered_map

///////////////////////////////////////////////////////////////////////////////
// Pleae note that the read_only_multimap offers a significant
// reduction in both run-time (approx * 0.4) and memory (approx * 0.6) 
// compared to std::unordered_multimap.
//#define MULTIMAP	std::unordered_multimap 
#define MULTIMAP	read_only_multimap

#define	DEFAULT_SCORE_WEIGHT	1.0f

//#define MAP	std::map
//#define MULTIMAP	std::multimap 

typedef enum {
	Seq_strand_unknown = 0,
	Seq_strand_plus = (1 << 0), 
	Seq_strand_minus = (1 << 1), 
	Seq_strand_both = Seq_strand_plus | Seq_strand_minus
	} Strand;

struct WordMatch
{
	unsigned int index;
	
	// The *effective* location of the 5' end of a word (offset to 
	// account for the fragment location *within* a word).
	int loc;
	Strand s;
		
	WordMatch()
	{
	};
	
	WordMatch(const unsigned int &m_index, const unsigned int &m_loc, 
		const Strand &m_strand):
		index(m_index), loc(m_loc), s(m_strand)
	{
	};
	
	// Return the 5' and 3' coordinates relative to the *plus* strand of 
	// the matching template sequence. A fixed coordinate system is needed
	// for computing binding site overlap and amplicon length.
	// The choice of plus target strand is arbitrary.
	int template_loc5(const int &m_start, const int &m_stop) const
	{
		if(s == Seq_strand_plus){
			return loc + m_start;
		}
		
		// Seq_strand_minus
		return loc - m_stop;
	};
	
	int template_loc3(const int &m_start, const int &m_stop) const
	{
		if(s == Seq_strand_plus){
			return loc + m_stop;
		}
		
		// Seq_strand_minus
		return loc - m_start;
	};
};

class Sequence
{	
	private:
		// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
		// ensure that structure variable are correctly serialized.
		#define SEQUENCE_MEMBERS \
			VARIABLE(std::string, def) \
			VARIABLE(std::deque<unsigned char>, seq_buffer) \
			VARIABLE(size_t, seq_len) \
			VARIABLE(float, score_weight) \
			VARIABLE(bool, is_active)
		
		#define VARIABLE(A, B) A B;
			SEQUENCE_MEMBERS
		#undef VARIABLE
	
		float extract_weight(const std::string &m_buffer) const;
	public:
	
		class const_iterator
		{
			private:
			
				std::deque<unsigned char>::const_iterator seq_iter;
				size_t seq_index;
			public:
				const_iterator(const std::deque<unsigned char>::const_iterator &m_seq_iter,
					const size_t &m_seq_index) :
					seq_iter(m_seq_iter), seq_index(m_seq_index)
				{
				};
				
				inline const_iterator& operator++()
				{
					++seq_index;
					
					if( !( seq_index & size_t(1) ) ){
						++seq_iter;
					}
					
					return *this;
				};
				
				inline char operator*() const
				{
					return ( ( seq_index & size_t(1) ) ? 
						(*seq_iter) : 
						(*seq_iter >> BITS_PER_BASE) ) & 0xF;
				};
				
				inline bool operator==(const const_iterator &m_rhs) const
				{
					return (seq_index == m_rhs.seq_index);
				};
				
				inline bool operator!=(const const_iterator &m_rhs) const
				{
					return (seq_index != m_rhs.seq_index);
				};
		};
		
		Sequence()
		{
			score_weight = DEFAULT_SCORE_WEIGHT;
			seq_len = 0;
			is_active = true;
		};
		
		Sequence(const std::string &m_seq, const float &m_weight = DEFAULT_SCORE_WEIGHT)
		{
			is_active = true;
			score_weight = m_weight;
			*this = m_seq;
		};
		
		Sequence(const std::string &m_def, const std::string &m_seq, const float &m_weight = DEFAULT_SCORE_WEIGHT)
		{
			is_active = true;
			score_weight = m_weight;
			def = m_def;
			*this = m_seq;
		};
		
		Sequence& operator=(const std::string &m_seq);
		Sequence& operator=(const std::deque<char> &m_seq);
		Sequence& operator+=(const std::deque<char> &m_seq);
		
		inline float weight() const
		{
			return score_weight;
		};
		
		inline void weight(const float &m_weight)
		{
			score_weight = m_weight;
		};
		
		inline bool active() const
		{
			return is_active;
		};
		
		inline void active(bool m_active)
		{
			is_active = m_active;
		};
		
		inline std::string defline() const
		{
			return def;
		};
		
		inline void defline(const std::string &m_def)
		{
			def = m_def;
			
			// Set the weight from the defline (if possible)
			score_weight = extract_weight(def);
			
			if(score_weight < 0.0f){
				throw __FILE__ ":Sequence::defline: Negative weights are not allowed!";
			}
		};
		
		// The actual sequence length (as opposed to the packed sequence length, since the 
		// data is stored as two bases per byte
		inline size_t length() const
		{
			return seq_len;
		};
		
		inline bool empty() const
		{
			return seq_buffer.empty();
		};
		
		Word subword(unsigned int m_loc, const unsigned int &m_len) const;
		
		// Since we pack words with less than or equal to m_cull_N_threshold N's, the default
		// value of WORD_LENGTH effectively disables culling.
		void pack(MULTIMAP<Word, WordMatch> &m_dbase, 
			const unsigned int &m_index, const unsigned int &m_degen_pack_threshold, 
			const float &m_min_gc, const float &m_max_gc,
			const unsigned int &m_min_oligo_length) const;
				
		inline unsigned char operator[](const unsigned int &m_index) const
		{
			const unsigned char b = seq_buffer[m_index/2];
			
			return (m_index%2 == 1) ? (b & 0xF): ( (b >> BITS_PER_BASE) & 0xF );
		};
		
		// Split a sequence by adding Base::EOS at the specified location
		inline void split_sequence(const unsigned int &m_index)
		{
			unsigned char b = seq_buffer[m_index/2];
			
			if(m_index%2 == 1){
				b = (b & ~0xF) | Base::EOS;
			}
			else{
				b = (b & ~0xF0) | (Base::EOS << BITS_PER_BASE);
			}

			seq_buffer[m_index/2] = b;
		};

		bool has_split(int m_loc, const int &m_len) const;

		const_iterator begin() const
		{
			return const_iterator(seq_buffer.begin(), size_t(0) );
		};
		
		const_iterator end() const
		{
			return const_iterator(seq_buffer.end(), seq_len);
		};
		
		// Append m_num_pad EOS symbols to the end of the sequence
		inline void pad(size_t m_num_pad)
		{
			while(m_num_pad > 0){
				
				if(seq_len%2 == 0){
					seq_buffer.push_back( (unsigned char)(0) );
				}
				
				++seq_len;
				--m_num_pad;
			}
		};
		
		friend std::ostream& operator << (std::ostream &m_s, const Sequence &m_seq);
		friend std::ostream& operator << (std::ostream &m_s, const Word &m_w);
		
		template<class T> friend size_t mpi_size(const T &m_obj);
		template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
			const T &m_obj);
		template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
			T &m_obj);
};

template<> size_t mpi_size<Sequence>(const Sequence &m_obj);
template<> unsigned char* mpi_pack<Sequence>(unsigned char* m_ptr, const Sequence &m_obj);
template<> unsigned char* mpi_unpack<Sequence>(unsigned char* m_ptr, Sequence &m_obj);

std::ostream& operator << (std::ostream &m_s, const Sequence &m_seq);

#endif // __SEQUENCE
