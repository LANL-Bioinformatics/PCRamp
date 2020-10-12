#include "assay.h"

#include <list>
#include <unordered_set>

using namespace std;

void select_words(MULTIMAP<Word, WordMatch> &m_dst, const MULTIMAP<Word, WordMatch> &m_src,
	const vector<PCR> &m_candidates, const bool &m_optimize_5, const bool &m_optimize_3,
	const float &m_threshold)
{

	if( m_src.empty() || m_candidates.empty() ){
		return;
	}
	
	// Search the source word keys (as opposed to the m_src multimap) to
	// avoid an inefficient search when the source data has lots of 
	// repeated words (i.e. eukaryotic genomes)
	const vector<Word> src_words = keys(m_src);
	
	// Extract all of the oligo words from the candidate assays
	deque< pair<Word /*oligo seq*/, unsigned int /*oligo size*/> > candidate_words;
					
	for(vector<PCR>::const_iterator i = m_candidates.begin();i != m_candidates.end();++i){
		
		for(Oligo *o_ptr = i->assay_oligos();*o_ptr != LAST_OLIGO;++o_ptr){
			
			// Defer the calculation of the oligo size (which is expensive)
			candidate_words.push_back( make_pair( i->oligo(*o_ptr), 0) );
							
			if(m_optimize_5 || m_optimize_3){
			
				// To enable growing and trimming during the assay search, we need to 
				// identify source words that match all of the allowed offsets of
				// the assay oligos.
				// ------XXXXXXX------- Original
				// -------XXXXXXX------ Shift right
				// --------XXXXXXX----- Shift right
				// ...
				// -----XXXXXXX-------- Shift left
				// ----XXXXXXX--------- Shift left
				// ...
				const Word &w = candidate_words.back().first;

				const int candidate_start = w.start();
				const int candidate_stop = w.stop();
				const int max_candidate_stop = w.max_size() - 1;

				if( m_optimize_5 && (candidate_start > 0) ){

					Word tmp(w);

					for(int j = 0;j < candidate_start;++j){

						tmp.shift_left();

						candidate_words.push_back( make_pair(tmp, 0) );
					}
				}

				if( m_optimize_3 && (candidate_stop < max_candidate_stop) ){

					Word tmp(w);

					for(int j = candidate_stop;j < max_candidate_stop;++j){

						tmp.shift_right();

						candidate_words.push_back( make_pair(tmp, 0) );
					}
				}
			}
		}
	}
	
	const size_t num_candidate_words = candidate_words.size();
	
	#pragma omp parallel for
	for(size_t i = 0;i < num_candidate_words;++i){
		
		// Since Word::size() is an expensive function, we compute it in parallel
		candidate_words[i].second = candidate_words[i].first.size()*m_threshold;
	}
	
	unordered_set<size_t> matched_words;
	
	#pragma omp parallel
	{
		deque<size_t> local;

		#pragma omp for
		for(size_t i = 0;i < num_candidate_words;++i){

			const pair<Word, unsigned int> &ref = candidate_words[i];

			// A local buffer that stores the best matching source oligos
			// to the current candidate word
			list<size_t> buffer;
			unsigned int best_match = ref.second;

			// Could store just best (or top) matches to each candiate word ...
			for(vector<Word>::const_iterator j = src_words.begin();j != src_words.end();++j){

				const unsigned int num_match = ref.first & *j;

				if(num_match >= best_match){

					if( best_match < num_match){
						buffer.clear();
					}

					best_match = num_match;

					buffer.push_back( j - src_words.begin() );						
				}
			}

			for(list<size_t>::const_iterator j = buffer.begin();j != buffer.end();++j){
				local.push_back(*j);
			}
		}

		// Merge the thread-specific local matches
		#pragma omp critical
		for(deque<size_t>::const_iterator i = local.begin();i != local.end();++i){			
			matched_words.insert(*i);
		}
	}
	
	for(unordered_set<size_t>::const_iterator i = matched_words.begin();i != matched_words.end();++i){

		typedef MULTIMAP<Word, WordMatch>::const_iterator I;

		const pair<I, I> range = m_src.equal_range(src_words[*i]);

		m_dst.insert(range.first, range.second);
	}
}
