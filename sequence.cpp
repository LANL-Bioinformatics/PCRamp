#include "sequence.h"
#include <ostream>
#include <iostream>
#include <algorithm>
#include <list>
#include <math.h>

#include <mpi.h>

using namespace std;

Sequence& Sequence::operator=(const string &m_seq)
{
	// Remove existing sequence
	seq_buffer.clear();
	
	// Make room to store this sequence. The current base packing
	// scheme requires 4 bits per base (== 2 bases per byte)
	seq_len = m_seq.size();
	
	const size_t packed_len = seq_len/2 + (seq_len%2 == 1);
		
	seq_buffer.resize(packed_len, 0); // <-- initialize to zero == no base
	
	string::const_iterator seq_iter = m_seq.begin();
	
	for(deque<unsigned char>::iterator i = seq_buffer.begin();i != seq_buffer.end();++i){
		
		(*i) = base_to_bits(*seq_iter) << BITS_PER_BASE;
		
		++seq_iter;
		
		if( seq_iter != m_seq.end() ){
		
			(*i) |= base_to_bits(*seq_iter);
			++seq_iter;
		}
	}
	
	return *this;
}

Sequence& Sequence::operator=(const deque<char> &m_seq)
{
	// Remove existing sequence
	seq_buffer.clear();
	
	// Make room to store this sequence. The current base packing
	// scheme requires 4 bits per base (== 2 bases per byte)
	seq_len = m_seq.size();
	
	const size_t packed_len = seq_len/2 + (seq_len%2 == 1);
		
	seq_buffer.resize(packed_len, 0); // <-- initialize to zero == no base
	
	deque<char>::const_iterator seq_iter = m_seq.begin();
	
	for(deque<unsigned char>::iterator i = seq_buffer.begin();i != seq_buffer.end();++i){
		
		(*i) = base_to_bits(*seq_iter) << BITS_PER_BASE;
		
		++seq_iter;
		
		if( seq_iter != m_seq.end() ){
		
			(*i) |= base_to_bits(*seq_iter);
			++seq_iter;
		}
	}
	
	return *this;
}

// Append m_seq to the current sequence
Sequence& Sequence::operator+=(const deque<char> &m_seq)
{
	for(deque<char>::const_iterator i = m_seq.begin();i != m_seq.end();++i){
		
		if(seq_len % 2 == 1){
			seq_buffer.back() |= base_to_bits(*i);
		}
		else{
			seq_buffer.push_back(base_to_bits(*i) << BITS_PER_BASE);
		}
		
		++seq_len;
	}
	
	return *this;
}

void Sequence::pack(MULTIMAP<Word, WordMatch> &m_dbase, const unsigned int &m_index, 
	const unsigned int &m_degen_pack_threshold, const float &m_min_gc, const float &m_max_gc,
	const unsigned int &m_min_oligo_length) const
{	
	Word w;
	size_t curr_word_size = 0;
	
	deque<unsigned char>::const_iterator iter = seq_buffer.begin();
	
	// Are we filtering by GC content?
	const bool gc_filter = (m_min_gc > 0.0f) || (m_max_gc < 1.0f);
	list<unsigned char> gc;
	unsigned int num_gc = 0;
	const unsigned int gc_mask = Base::G | Base::C;
	const float norm = 1.0f/w.max_size();
	
	int loc;
	
	for(loc = 1;iter != seq_buffer.end();++loc){

		unsigned char b = 0;
		
		if(loc % 2 == 1){ // This formula depends on the starting value of loc
			b = (*iter >> BITS_PER_BASE) & 0xF;
		}
		else{
			b = *iter & 0xF;
			++iter;
		}
		
		w.push_back(b);
		
		// Only count "real" bases when computing the word size
		curr_word_size += (b != Base::EOS);
		
		if(gc_filter){
			
			if( gc.size() == w.max_size() ){
				
				num_gc -= ( (gc.front() & gc_mask) != 0 );
				
				gc.pop_front();
			}
			
			gc.push_back(b);
			
			num_gc += ( (b & gc_mask) != 0 );
			
			const float fraction_gc = num_gc*norm;
			
			if( (fraction_gc < m_min_gc) || (fraction_gc > m_max_gc) ){
				
				curr_word_size = min(curr_word_size, w.max_size() - 1UL);
				continue;
			}
		}
		
		if(w.degeneracy() > m_degen_pack_threshold){
		
			curr_word_size = min(curr_word_size, w.max_size() - 1UL);
			continue;
		}
		
		if( curr_word_size < w.max_size() ){
			
			if(curr_word_size >= m_min_oligo_length){
			
				Word tmp = w;

				tmp.center();

				// Save the location of the 5' end correctly offset for the
				// start of the word (which does not occupy the full max_size).
				m_dbase.insert( 
					make_pair( tmp, 
						WordMatch(m_index, loc - int(curr_word_size) - tmp.start(), 
							Seq_strand_plus) ) );

				// Since taking the complement "uncenters" a word, we need to recenter it
				// after computing the complement.
				tmp = tmp.complement();
				tmp.center();

				m_dbase.insert( 
					make_pair( tmp, 
						WordMatch(m_index, loc - 1 + tmp.start(), 
							Seq_strand_minus) ) );
			}
			
		}
		else{ //curr_word_size == w.max_size()
		
			// Save the location of the 5' end
			m_dbase.insert( 
				make_pair( w, 
					WordMatch(m_index, loc - curr_word_size, Seq_strand_plus) ) );
			
			m_dbase.insert( 
				make_pair( w.complement(), 
					WordMatch(m_index, loc - 1, Seq_strand_minus) ) );
					
			--curr_word_size;
		}
	}
		
	// Keep shifting the last word to the left so we capture *all* of the word offsets
	while(curr_word_size > 0){
		
		w.shift_left();
		
		--curr_word_size;
		
		if(gc_filter){
			
			if( gc.size() == w.max_size() ){
				
				num_gc -= ( (gc.front() & gc_mask) != 0 );
				
				gc.pop_front();
			}
			
			const float fraction_gc = num_gc*norm;
			
			if( (fraction_gc < m_min_gc) || (fraction_gc > m_max_gc) ){
				
				// Don't increment loc, since the 3' end of the word is
				// no longer changing (by shifting the word left, we are
				// change the location of the of 5' end of the word -- which
				// is not stored!).
				//++loc;
				continue;
			}
		}
		
		if(w.degeneracy() > m_degen_pack_threshold){
			
			// Don't increment loc, since the 3' end of the word is
			// no longer changing (by shifting the word left, we are
			// change the location of the of 5' end of the word -- which
			// is not stored!).
			//++loc;
			continue;
		}
		
		if(curr_word_size >= m_min_oligo_length){
		
			Word tmp = w;

			// Center a copy of the current word
			tmp.center();

			m_dbase.insert( 
				make_pair( tmp, 
					WordMatch(m_index, loc - 1 - int(curr_word_size) - tmp.start(), Seq_strand_plus) ) );
			
			// Since taking the complement "uncenters" a word, we need to recenter it
			// after taking the complement.
			tmp = tmp.complement();
			tmp.center();

			m_dbase.insert( 
				make_pair( tmp, 
					WordMatch(m_index, loc - 2 + tmp.start(), Seq_strand_minus) ) );
					
		}
		
		// Don't increment loc, since the 3' end of the word is
		// no longer changing (by shifting the word left, we are
		// change the location of the of 5' end of the word -- which
		// is not stored!).
		//++loc;
	}
	
	// Make the multimap iterable
	m_dbase.sort();
}

Word Sequence::subword(unsigned int m_loc, const unsigned int &m_len) const
{
	Word ret;
	
	if( (m_loc + m_len) > length() ){
		throw __FILE__ ":Sequence::subword: word is out of bounds";
	}
	
	deque<unsigned char>::const_iterator seq = seq_buffer.begin() + m_loc/2;
	
	for(unsigned int i = 0;i < m_len;++i, ++m_loc){
		
		unsigned char b = 0;
		
		if(m_loc % 2 == 0){
			b = (*seq >> BITS_PER_BASE) & 0xF;
		}
		else{
			b = *seq & 0xF;
			++seq;
		}
		
		// Multi-fasta record sequences can contain EOS symbols scattered throughout
		// the sequence (as padding between records). As a result, it is no longer 
		// an error to encounter an unexpected Base::EOS.
		//if(b == Base::EOS){
		//	throw __FILE__ ":Sequence::subword: Encountered EOS";
		//}
		
		ret.push_back(b);
	}
	
	return ret;
}

bool Sequence::has_split(int m_loc, const int &m_len) const
{
	if( ( (m_loc + m_len) > length() ) || (m_loc < 0) || (m_len < 0) ){
		throw __FILE__ ":Sequence::has_split: range is out of bounds";
	}
	
	deque<unsigned char>::const_iterator seq = seq_buffer.begin() + m_loc/2;
	
	for(unsigned int i = 0;i < m_len;++i, ++m_loc){
		
		unsigned char b = 0;
		
		if(m_loc % 2 == 0){
			b = (*seq >> BITS_PER_BASE) & 0xF;
		}
		else{
			b = *seq & 0xF;
			++seq;
		}
		
		if(b == Base::EOS){
			return true;
		}
	}
	
	return false;
}

float Sequence::extract_weight(const std::string &m_buffer) const
{
	// Extract user-defined weights in the form of:
	// 	[w=xxxx]
	// where xxxx is a floating point weight. Use a finite state
	// machine for matching.
	enum {NO_MATCH, LEFT, WEIGHT, EQUAL, VALUE, RIGHT};
	
	const string::size_type len = m_buffer.size();
	int match = NO_MATCH;
	
	string::size_type start = string::npos;
	string::size_type stop = string::npos;
	
	// This implementation is way to convoluted! I should just use the C++11 regex library ...
	for(string::size_type i = 0;i < len;++i){
		
		switch(match){
			case NO_MATCH:
				
				if(m_buffer[i] == '['){
					match = LEFT;
				}
				
				break;
			case LEFT:
				
				switch(m_buffer[i]){
					case 'W':
					case 'w':
						match = WEIGHT;
						break;
					case ' ':
					case '\t':
					case '[':
						// Skip spaces and multiple '['
						break;
					default:
						match = NO_MATCH;
						break;
				};
								
				break;
			case WEIGHT:
				
				switch(m_buffer[i]){
					case '=':
						match = EQUAL;
						break;
					case ' ':
					case '\t':
						// Skip spaces
						break;
					case '[':
						match = LEFT;
						break;
					default:
						match = NO_MATCH;
						break;
				};
				
				break;
			case EQUAL:
				
				switch(m_buffer[i]){
					case ' ':
					case '\t':
						// Skip spaces
						break;
					case '-':
					case '+':
					case '.':
					case 'e':
					case '0':
					case '1':
					case '2':
					case '3':
					case '4':
					case '5':
					case '6':
					case '7':
					case '8':
					case '9':
						start = i;
						match = VALUE;
						break;
					case '[':
						match = LEFT;
						break;
					default:
						match = NO_MATCH;
						break;
				};
				
				break;
			case VALUE:
				
				switch(m_buffer[i]){
					case ' ':
					case '\t':
					
						match = RIGHT;
						break;
					case ']':
					
						assert(stop >= start);	
						return atof( m_buffer.substr(start, stop - start + 1).c_str() );
					case '-':
					case '+':
					case '.':
					case 'e':
					case '0':
					case '1':
					case '2':
					case '3':
					case '4':
					case '5':
					case '6':
					case '7':
					case '8':
					case '9':
					
						stop = i;
						match = VALUE;
						break;
					case '[':
						match = LEFT;
						break;
					default:
						match = NO_MATCH;
						break;
				};
				
				break;
			case  RIGHT:
				
				switch(m_buffer[i]){
					case ' ':
					case '\t':
						// Skip spaces
						break;
					case ']':
					
						assert(stop >= start);
						return atof( m_buffer.substr(start, stop - start + 1).c_str() );
					case '[':
						match = LEFT;
						break;
					default:
						match = NO_MATCH;
						break;
				};
				
				break;
			default:
				throw __FILE__ ":Sequence::extract_weight: Unknown state";
		};
	}
	
	return DEFAULT_SCORE_WEIGHT;
	
}

ostream& operator << (ostream &m_s, const Sequence &m_seq)
{
	size_t index = 0;
	
	for(deque<unsigned char>::const_iterator i = m_seq.seq_buffer.begin();i != m_seq.seq_buffer.end();++i){
		
		// High order bits
		unsigned char b = (*i >> BITS_PER_BASE) & 0xF;
		
		m_s << bits_to_base(b);
		
		++index;
		
		// Low order bits
		b = *i & 0xF;
		
		if( index != m_seq.length() ){
			m_s << bits_to_base(b);
		}
		
		++index;
	}
	
	return m_s;
}
