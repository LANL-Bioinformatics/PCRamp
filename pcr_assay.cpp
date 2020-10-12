#include <math.h>
#include <algorithm>
#include <unordered_set>
//#include <set>
#include "assay.h"

#define SET	std::unordered_set
//#define SET	std::set

using namespace std;

void PCR::collect_candidates(deque<PCROligos> &m_amplicons, 
	MAP<unsigned int, float> &m_f_identity, 
	MAP<unsigned int, float> &m_r_identity,
	const vector<Word> &m_keys, 
	const MULTIMAP<Word, WordMatch> &m_db, const deque<Sequence> &m_seq,
	const float &m_threshold,
	const std::pair<int, int> &m_amplicon_range) const
{	
	m_amplicons.clear();
		
	// Step 1: Extract the sequences targeted by the input assay
	// Start by finding the matches to the forward and reverse oligos.
	// Rather than storing words, store indicies into the key array
	deque<unsigned int> f_match;
	deque<unsigned int> r_match;
	
	// The PCR == sqrt(f*r) >= m_threshold
	// If f == 1 (best case), the smallest allowed value of r = m_threshold^2.
	// Same argument leads to the smallest allowed value of f = m_threshold^2.
	match_words( f_match, oligo(FORWARD), m_keys, m_threshold*m_threshold );	
	match_words( r_match, oligo(REVERSE), m_keys, m_threshold*m_threshold );
		
	///////////////////////////////////////////////////////////////////////////////////
	// Find the set of sequences matched by {F(+), R(-)}
	///////////////////////////////////////////////////////////////////////////////////
	deque<OligoMatch> o_match;

	find_oligo_match(o_match, f_match, FORWARD, Seq_strand_plus, m_keys, m_db, m_seq);
	find_oligo_match(o_match, r_match, REVERSE, Seq_strand_minus, m_keys, m_db, m_seq);
		
	// Sort by sequence id, then sort by oligo location
	sort( o_match.begin(), o_match.end() );
		
	find_amplicon_match(m_amplicons, o_match, FORWARD, REVERSE, m_seq, m_amplicon_range);
	
	///////////////////////////////////////////////////////////////////////////////////
	// Find the set of sequences matched by {R(+), F(-)}
	///////////////////////////////////////////////////////////////////////////////////
	o_match.clear();

	find_oligo_match(o_match, f_match, FORWARD, Seq_strand_minus, m_keys, m_db, m_seq);
	find_oligo_match(o_match, r_match, REVERSE, Seq_strand_plus, m_keys, m_db, m_seq);

	// First sort by sequence id, then sort by oligo location
	sort( o_match.begin(), o_match.end() );
		
	find_amplicon_match(m_amplicons, o_match, REVERSE, FORWARD, m_seq, m_amplicon_range);
	
	m_f_identity.clear();
	m_r_identity.clear();

	// Find the keys that we will need to evaluate
	for(deque<PCROligos>::const_iterator i = m_amplicons.begin();i != m_amplicons.end();++i){

		m_f_identity[i->f] = 0.0f;
		m_r_identity[i->r] = 0.0f;
	}
}

void PCR::collect_multiplex_background_candidates(
	const std::vector<Word> &m_keys, 
	const MULTIMAP<Word, WordMatch> &m_db, const std::deque<Sequence> &m_seq, const Options &m_opt)
{
	if( m_keys.empty() ){
		return;
	}
	
	// Check for *single* primer overlap with a multiplex background amplicon
	multiplex_background_f_identity.clear();
	multiplex_background_r_identity.clear();
	
	deque<unsigned int> primer_match;
	
	// Search each oligo independently
	match_words(primer_match, oligo(FORWARD), m_keys, m_opt.background_threshold);	
	
	// Find the keys that we will need to evaluate
	for(deque<unsigned int>::const_iterator i = primer_match.begin();i != primer_match.end();++i){
		multiplex_background_f_identity[*i] = 0.0f;
	}
	
	primer_match.clear();
	
	// Search each oligo independently
	match_words(primer_match, oligo(REVERSE), m_keys, m_opt.background_threshold);	
	
	// Find the keys that we will need to evaluate
	for(deque<unsigned int>::const_iterator i = primer_match.begin();i != primer_match.end();++i){
		multiplex_background_r_identity[*i] = 0.0f;
	}
}
			
void  PCR::sort_candidates(std::deque<PCROligos> &m_amplicons, 
	const MAP<unsigned int, float> &m_f_identity, 
	const MAP<unsigned int, float> &m_r_identity)
{
	sort( m_amplicons.begin(), m_amplicons.end(), sort_by_PCR_score(m_f_identity, m_r_identity) );
}

bool PCR::maximize_degeneracy(const vector<Word> &m_keys, NucCruc &m_melt, const Options &m_opt)
{
	for(deque<PCROligos>::const_iterator i = target_amplicons.begin();i != target_amplicons.end();++i){
		
		const Word local_f = f | m_keys[i->f];
		const Word local_r = r | m_keys[i->r];
	
		// Select degenerate oligos that are individually valid (i.e. have valid degeneracies
		// and hairpin melting temperatures)
		if( (local_f.degeneracy() <= m_opt.degen) && is_valid(FORWARD, local_f, m_melt, m_opt, true /*check dimer*/) ){
			f = local_f;
		}
		
		if( (local_r.degeneracy() <= m_opt.degen) && is_valid(REVERSE, local_r, m_melt, m_opt, true /*check dimer*/) ){
			r = local_r;
		}
	}

	// If needed, *decrease* the degeneracy of some (or all) of the assay oligos to make sure that
	// there is no illegal heterodimer formation.
	float min_dimer_tm = max_dimer_tm(m_melt, m_opt);
	
	const int last_f = f.stop();
	const int last_r = r.stop();
	
	// Find the single base reduction in degeneracy that yeilds the largest decrease in the max dimer melting temperature
	while(min_dimer_tm > m_opt.max_dimer){
		
		// Start with An impossibly large melting temperature so we always accept a lower degeneracy
		// assay
		float curr_dimer_tm = 1.0e6f;
		Oligo best_oligo = LAST_OLIGO;
		Word best;
		
		// Search all possible forward primer reductions in degeneracy
		for(int i = f.start();i <= last_f;++i){
		
			const unsigned char curr = f.get(i);

			for(unsigned char b = Base::A;b <= Base::T;b <<= 1){

				// Skip bases that do not include "b" in combination with another
				// base
				const unsigned char d = (curr & ~b);

				// if d == 0, we have masked out a single base (no degeneracy)
				// if d == curr, we have masked a base that does not appear in curr
				if( !d || (d == curr) ){
					continue;
				}

				// Decrease the degeneracy by removing b
				f.unmask(b, i);

				const float tm = max_dimer_tm(m_melt, m_opt);
				
				if(tm < curr_dimer_tm){
					
					curr_dimer_tm = tm;
					best_oligo = FORWARD;
					best = f;
				}
				
				// Restore the original base
				f.mask(b, i);
			}			
		}
		
		// Search all possible reverse primer reductions in degeneracy
		for(int i = r.start();i <= last_r;++i){
		
			const unsigned char curr = r.get(i);

			for(unsigned char b = Base::A;b <= Base::T;b <<= 1){

				// Skip bases that do not include "b" in combination with another
				// base
				const unsigned char d = (curr & ~b);

				// if d == 0, we have masked out a single base (no degeneracy)
				// if d == curr, we have masked a base that does not appear in curr
				if( !d || (d == curr) ){
					continue;
				}

				// Decrease the degeneracy by removing b
				r.unmask(b, i);

				const float tm = max_dimer_tm(m_melt, m_opt);
				
				if(tm < curr_dimer_tm){
					
					curr_dimer_tm = tm;
					best_oligo = REVERSE;
					best = r;
				}
				
				// Restore the original base
				r.mask(b, i);
			}			
		}
		
		switch(best_oligo){
			case FORWARD:
				f = best;
				break;
			case REVERSE:
				r = best;
				break;
			default:
				// The greedy reduction in oligo degeneracy has lead to a set of sequences that
				// are non-degenerate, but do not satisfy the melting temperature constraints.
				return false;
		};
		
		min_dimer_tm = curr_dimer_tm;
	}
	
	return true;
}

float PCR::max_dimer_tm(NucCruc &m_melt, const Options &m_opt) const
{
	float ret = 0.0f;
	
	const double degen_f = oligo(FORWARD).degeneracy();
	const double degen_r = oligo(REVERSE).degeneracy();
	
	/////////////////////////////////////////////////////////////////////////////////////////
	// Check for heterodimer formation
	/////////////////////////////////////////////////////////////////////////////////////////
	Word f_iter = oligo(FORWARD).begin();

	m_melt.strand(m_opt.primer_strand/degen_f, m_opt.primer_strand/degen_r);

	do{ // Enumerate possibly degenerate oligos

		const string local_f = f_iter.str();

		m_melt.set_query(local_f);

		Word r_iter = oligo(REVERSE).begin();

		do{ // Enumerate possibly degenerate oligos

			const string local_r = r_iter.str();

			m_melt.set_target(local_r);

			const float tm = m_melt.approximate_tm_heterodimer();

			ret = max(ret, tm);

		} while( oligo(REVERSE).next(r_iter) );

	} while( oligo(FORWARD).next(f_iter) );

	return ret;
}

float PCR::compute_coverage(const deque<PCROligos> &m_amplicons, 
	const MAP<unsigned int, float> &m_f, 
	const MAP<unsigned int, float> &m_r, 
	const float &m_threshold) const
{
	if( m_amplicons.empty() ){
		return 0;
	}
	
	// Compute the score by summing the weights of each detected sequence
	double ret = 0.0; // Accumulate in double
	
	SET<unsigned int> valid;
		
	for(deque<PCROligos>::const_iterator i = m_amplicons.begin();i != m_amplicons.end();++i){

		MAP<unsigned int, float>::const_iterator f_iter = m_f.find(i->f);
		MAP<unsigned int, float>::const_iterator r_iter = m_r.find(i->r);

		assert( ( f_iter != m_f.end() ) && ( r_iter != m_r.end() ) );

		const float local = sqrtf(f_iter->second * r_iter->second);

		if( (local >= m_threshold) && ( valid.find(i->index) == valid.end() ) ){
		
			valid.insert(i->index);
			ret += i->weight;
		}
	}
	
	return ret;
}

float PCR::compute_multiplex_background_coverage(const float &m_threshold) const
{
	if( multiplex_background_f_identity.empty() && multiplex_background_r_identity.empty() ){
		return 0;
	}
	
	// Compute the score by summing the weights of each detected sequence
	double ret = 0.0; // Accumulate in double
	
	SET<unsigned int> valid;
		
	for(MAP<unsigned int, float>::const_iterator i = multiplex_background_f_identity.begin();
		i != multiplex_background_f_identity.end();++i){

		if( (i->second >= m_threshold) && ( valid.find(i->first) == valid.end() ) ){
		
			valid.insert(i->first);
			ret += 1.0; // Each multiplex background gets weight 1.0
		}
	}
	
	for(MAP<unsigned int, float>::const_iterator i = multiplex_background_r_identity.begin();
		i != multiplex_background_r_identity.end();++i){

		if( (i->second >= m_threshold) && ( valid.find(i->first) == valid.end() ) ){
		
			valid.insert(i->first);
			ret += 1.0; // Each multiplex background gets weight 1.0
		}
	}
	
	return ret;
}

void PCR::find_amplicon_match(deque<PCROligos> &m_amplicons, 
	deque<OligoMatch> &m_match, const Oligo &m_plus_oligo, const Oligo &m_minus_oligo, 
	const deque<Sequence> &m_seq, const pair<int, int> &m_amplicon_range) const
{
	const int plus_start = oligo(m_plus_oligo).start();
	const int plus_stop = oligo(m_plus_oligo).stop();
	
	const int minus_start = oligo(m_minus_oligo).start();
	const int minus_stop = oligo(m_minus_oligo).stop();
		
	for(deque<OligoMatch>::const_iterator plus = m_match.begin();plus != m_match.end();++plus){
		
		if(plus->o != m_plus_oligo){
			continue;
		}
				
		for(deque<OligoMatch>::const_iterator minus = plus;minus != m_match.end();++minus){

			// Only search minus primer binding sites that target the same sequence
			// as the plus primer
			if(plus->index != minus->index){
				break;
			}

			if(minus->o != m_minus_oligo){
				continue;
			}

			// Skip minus primers that overlap the plus primer
			if( plus->template_loc3(plus_start, plus_stop) 
				>= minus->template_loc5(minus_start, minus_stop) ){
				continue;
			}

			int amp_start = plus->template_loc5(plus_start, plus_stop);

			// If mismatches are allowed, the primer could dangle off the end
			// of the template sequence.
			const int amp_stop = min( int( minus->template_loc3(minus_start, minus_stop) ),
				int(m_seq[plus->index].length() - 1) );

			int amp_len = amp_stop - amp_start + 1;

			// Skip reverse primers that generate an amplicon that is too small
			if(amp_len < m_amplicon_range.first){
				continue;
			}

			// Stop searching reverse primers when the amplicon length is too large,
			// such that we run off the end of the sequence or exceed the 
			// maximum allowed amplicon length
			if(amp_len > m_amplicon_range.second){
				break;	
			}

			//bool splitter = false;

			//try{
			//	splitter = m_seq[plus->index].has_split(amp_start, amp_len);
			//}
			//catch(...){
				// DEBUG
			//	cerr << "seq index = " << plus->index << endl;
			//	cerr << "defline = " << m_seq[plus->index].defline() << endl;
			//	cerr << "\tamp_start = " << amp_start << endl;
			//	cerr << "\tamp_len = " << amp_len << endl;
			//	cerr << "\tplus oligo = " << oligo(m_plus_oligo).str() << endl;
			//	cerr << "\tminus oligo = " << oligo(m_minus_oligo).str() << endl;
			//}

			// If we allow mismatches, there is the possiblity that a primer can dangle
			// off the 5' of a sequence
			if(amp_start < 0){

				amp_len += amp_start;
				amp_start = 0;
			}

			if( m_seq[plus->index].has_split(amp_start, amp_len) ){
				break;
			}

			//if( m_seq[plus->index].has_split(amp_start, amp_len) ){
			//	break;
			//}

			// This is a valid candidate amplicon!
			if(plus->o == FORWARD){
				m_amplicons.push_back(
					PCROligos( (unsigned int)(plus->index), 
						m_seq[plus->index].weight(),
						plus->key_index, 
						minus->key_index) );
			}
			else{
				m_amplicons.push_back(
					PCROligos( (unsigned int)(plus->index), 
						m_seq[plus->index].weight(),
						minus->key_index, 
						plus->key_index) );
			}
		}
	}
}

void PCR::extract_amplicon_seq(deque<string> &m_amplicons, 
	deque<OligoMatch> &m_match, const Oligo &m_plus_oligo, const Oligo &m_minus_oligo, 
	const deque<Sequence> &m_seq, const pair<int, int> &m_amplicon_range,
	deque<AmpliconBounds> *m_bounds_ptr /*= NULL*/)
{
	const int plus_start = oligo(m_plus_oligo).start();
	const int plus_stop = oligo(m_plus_oligo).stop();
	
	const int minus_start = oligo(m_minus_oligo).start();
	const int minus_stop = oligo(m_minus_oligo).stop();
	
	for(deque<OligoMatch>::const_iterator plus = m_match.begin();plus != m_match.end();++plus){
		
		if(plus->o != m_plus_oligo){
			continue;
		}
	
		for(deque<OligoMatch>::const_iterator minus = plus;minus != m_match.end();++minus){

			// Only search minus primer binding sites that target the same sequence
			// as the plus primer
			if(plus->index != minus->index){
				break;
			}

			if(minus->o != m_minus_oligo){
				continue;
			}
		
			// Skip minus primers that overlap the plus primer
			if( plus->template_loc3(plus_start, plus_stop) 
				>= minus->template_loc5(minus_start, minus_stop) ){
				continue;
			}

			const int amp_len = minus->template_loc3(minus_start, minus_stop) - 
				plus->template_loc5(plus_start, plus_stop) + 1;

			// Skip reverse primers that generate an amplicon that is too small
			if(amp_len < m_amplicon_range.first){
				continue;
			}

			// Stop searching reverse primers when the amplicon length is too large
			if(amp_len > m_amplicon_range.second){
				break;
			}
			
			// The non-primer amplicon length skips the primer binding sites (since these
			// sites are valid targets for later assays). Include MULTIPLEX_AMPLICON_PADDING 
			// extra bases on both the 5' and 3' ends of the amplicon sequence.
			const unsigned int amp_start = 
				plus->template_loc3(plus_start, plus_stop) + 1 - MULTIPLEX_AMPLICON_PADDING;
			
			const unsigned int non_primer_amp_len = 
				minus->template_loc5(minus_start, minus_stop) - 
				amp_start + 2*MULTIPLEX_AMPLICON_PADDING;

			string amp_seq(non_primer_amp_len, '?');
			
			const Sequence &template_seq = m_seq[plus->index];
			bool valid_amp = true;

			for(unsigned int i = 0;i < non_primer_amp_len;++i){

				const unsigned char b = template_seq[amp_start + i];

				// Make sure any potential amplicon does not span multiple sequences
				if(b == Base::EOS){

					valid_amp = false;
					break;
				}

				amp_seq[i] = bits_to_base(b);
			}
			
			if(!valid_amp){

				// This primer pair spans two or more sequences!
				break;
			}

			m_amplicons.push_back(amp_seq);

			if(m_bounds_ptr != NULL){

				// Store the location of the begining and end of amplicon on the target sequence
				// (including the forward and reverse primers).
				// For successfull assays, we will edit the target sequences by adding
				// Base::IOS symbols at the first and last base of the amplicon
				// to effectively "break" the sequence and prevent future PCR primer
				// pairs from flanking this current set of primer pairs.
				m_bounds_ptr->push_back( AmpliconBounds(plus->index, 
					plus->template_loc5(plus_start, plus_stop), 		// Begin
					minus->template_loc3(minus_start, minus_stop) ) );  // End
			}
		}
	}
}

void PCR::find_target_match(BitSet &m_match, const vector<Word> &m_keys, 
	const MULTIMAP<Word, WordMatch> &m_db, const deque<Sequence> &m_seq, const Options &m_opt)
{	
	// Don't clear assay matches -- append as needed
	//m_match.clear();
	
	if( m_match.size() != m_seq.size() ){
		m_match.resize( m_seq.size(), false );
	}
	
	collect_candidates(target_amplicons, target_f_identity, target_r_identity,
		m_keys, m_db, m_seq, 
		m_opt.target_threshold,
		m_opt.target_amplicon_range);
	
	if( target_amplicons.empty() ){
		return;
	}
		
	update_target_candidates(m_keys, m_opt.use_taq_mama);
		
	for(deque<PCROligos>::const_iterator i = target_amplicons.begin();i != target_amplicons.end();++i){

		MAP<unsigned int, float>::const_iterator f_iter = target_f_identity.find(i->f);
		MAP<unsigned int, float>::const_iterator r_iter = target_r_identity.find(i->r);

		assert( ( f_iter != target_f_identity.end() ) && ( r_iter != target_r_identity.end() ) );

		const float local = sqrtf(f_iter->second * r_iter->second);

		if(local >= m_opt.target_threshold){
			m_match[i->index] = true;
		}
	}
}

void PCR::random_assay(const deque<Sequence> &m_seq, 
	NucCruc &m_melt, const Options &m_opt, unsigned int &m_seed, ostream &m_out)
{
	const unsigned int max_sequence_iter = 100;
	const unsigned int max_assay_iter = 100;
	
	assert( !m_seq.empty() );

	// Keep iterating until we find a sequence that yeilds a valid assay	
	unsigned int sequence_iteration = 0;
	
	deque<size_t> indicies;
	const size_t num_seq = m_seq.size();
	
	for(size_t i = 0;i < num_seq;++i){
	
		if( m_seq[i].active() ){
			indicies.push_back(i);
		}
	}
	
	const size_t num_active = indicies.size();
	
	if(num_active == 0){
		throw __FILE__ ":PCR::random_assay: No active sequences found";
	}
	
	while(true){
		
		++sequence_iteration;
		
		if(sequence_iteration > max_sequence_iter){
			throw __FILE__ ":PCR::random_assay: Unable to generate a valid initial assay to test!";
		}

		//const Sequence& target_seq = m_seq[ indicies[random()%indicies.size()] ];
		//const Sequence& target_seq = m_seq[ random_index(m_cumulative_weight, m_seed) ];

		const Sequence& target_seq = m_seq[ indicies[rand_r(&m_seed)%num_active] ];
		
		const int len = target_seq.length();

		if(len < m_opt.target_amplicon_range.first){
			throw __FILE__ ":PCR::random_assay: sequence length is too small!";
		}

		unsigned int assay_iteration = 0;
		
		while(true){
			
			++assay_iteration;

			if(assay_iteration > max_assay_iter){
				break;
			}

			// Randomly select the oligo length parameters:
			const int f_len = m_opt.primer_range.first + rand_r(&m_seed)%(m_opt.primer_range.second - m_opt.primer_range.first + 1);
			const int r_len = m_opt.primer_range.first + rand_r(&m_seed)%(m_opt.primer_range.second - m_opt.primer_range.first + 1);
			
			if( (f_len + r_len) > len){
				continue;
			}
			
			const int f_start = random_location(0, (len + 1) - m_opt.target_amplicon_range.first, &m_seed);
			
			oligo( FORWARD, target_seq.subword(f_start, f_len) );
			
			// Sequence::subword can return an incomplete oligo if this is a padded, multi-record sequence
			if( oligo(FORWARD).size() != size_t(f_len) ){
				continue;
			}
			
			const double degen_f = oligo(FORWARD).degeneracy();
			
			if(degen_f > m_opt.degen){
				continue;
			}
			
			/////////////////////////////////////////////////////////////////////////////////////////
			// Check for perfect match Tm as well as hairpin and homodimer formation
			/////////////////////////////////////////////////////////////////////////////////////////
			if( !is_valid(FORWARD, oligo(FORWARD), m_melt, m_opt, true /* check dimer */) ){
				continue;
			}
			
			// The following line of code introduced a very subtle bug. By randomly selecting the
			// start of the reverse primer from a large range (from the end of the forward primer to the end
			// of the source sequence) we preferentially selected for assay that were near the 3' end of a 
			// sequence (since a larger fraction of these assays would have a valid amplicon length). Instead,
			// restrict the random reverse primer location to be from the end of the forward primer to the
			// maximum amplicon length.
			//const int r_start = random_location(coverage, f_start + f_len, len - r_len); <-- don't use!!!
			//const int r_start = random_location(f_start + f_len, 
			//	min( (len + 1) - r_len, (f_start + m_opt.target_amplicon_range.second + 1) - r_len),
			//	&m_seed);
			
			const int r_start = random_location(f_start + m_opt.target_amplicon_range.first - r_len, 
				min( (len + 1) - r_len, (f_start + m_opt.target_amplicon_range.second + 1) - r_len),
				&m_seed);
			
			const int amp_len = r_start - f_start + r_len;
			
			if( (amp_len > m_opt.target_amplicon_range.second) || 
			    (amp_len < m_opt.target_amplicon_range.first) ){
				continue;
			}
			
			oligo( REVERSE, target_seq.subword(r_start, r_len).complement() );
			
			// Sequence::subword can return an incomplete oligo if this is a padded, multi-record sequence
			if( oligo(REVERSE).size() != size_t(r_len) ){
				continue;
			}
			
			const double degen_r = oligo(REVERSE).degeneracy();
			
			if(degen_r > m_opt.degen){
				continue;
			}
			
			// Make sure that this primer pair does not span two or more distinct sequences
			// within the length of the amplicon
			if( target_seq.has_split(f_start, amp_len) ){
				continue;
			}

			/////////////////////////////////////////////////////////////////////////////////////////
			// Check for perfect match Tm as well as hairpin and homodimer formation
			/////////////////////////////////////////////////////////////////////////////////////////
			if( !is_valid(REVERSE, oligo(REVERSE), m_melt, m_opt, true /* check dimer */) ){
				continue;
			}
			
			// Check for heterodimer formation
			if(max_dimer_tm(m_melt, m_opt) > m_opt.max_dimer){
				continue;
			}
			
			center();
			
			if(m_opt.output_filter > Options::VERBOSE){
				
				// If we get here, we have a valid assay
				m_out << "Init assay (tried " << sequence_iteration << " seq and " << assay_iteration << " assays): ";
			
				write(m_out);
			
				m_out << " : Amplicon length = " << amp_len << endl;
			}
			
			return;
		}
	}	
}

float PCR::compute_oligo_overlap(const deque<PCR> &m_pool) const
{
	float best_f = 0.0;
	float best_r = 0.0;
	
	for(deque<PCR>::const_iterator i = m_pool.begin();i != m_pool.end();++i){
		
		// Since forward and reverse primers are interchangable, we need to 
		// compute the similarity to both primers
		best_f = max( best_f, oligo(FORWARD).max_overlap( i->oligo(FORWARD) ) );
		best_f = max( best_f, oligo(FORWARD).max_overlap( i->oligo(REVERSE) ) );
		
		best_r = max( best_r, oligo(REVERSE).max_overlap( i->oligo(FORWARD) ) );
		best_r = max( best_r, oligo(REVERSE).max_overlap( i->oligo(REVERSE) ) );
	}
	
	return ( (best_f == 1.0f) ? MULTIPLEX_OLIGO_REUSE_BONUS : best_f ) + 
	       ( (best_r == 1.0f) ? MULTIPLEX_OLIGO_REUSE_BONUS : best_r );
}

deque<Sequence> PCR::collect_unique_amplicons(
	const vector<Word> &m_keys, 
	const MULTIMAP<Word, WordMatch> &m_db, 
	const deque<Sequence> &m_seq, 
	const float &m_threshold,
	const pair<int, int> &m_amplicon_range,
	std::deque<AmpliconBounds> *m_bounds_ptr /*=NULL*/)
{	
	deque<string> amplicons;
	
	// Step 1: Extract the sequences targeted by the input assay
	// Start by finding the matches to the forward and reverse oligos.
	// Rather than storing words, store indicies into the key array
	deque<unsigned int> f_match;
	deque<unsigned int> r_match;
	
	// The PCR == sqrt(f*r) >= m_threshold
	// If f == 1 (best case), the smallest allowed value of r = m_threshold^2.
	// Same argument leads to the smallest allowed value of f = m_threshold^2.
	match_words( f_match, oligo(FORWARD), m_keys, m_threshold*m_threshold );	
	match_words( r_match, oligo(REVERSE), m_keys, m_threshold*m_threshold );
		
	///////////////////////////////////////////////////////////////////////////////////
	// Find the set of sequences matched by {F(+), R(-)}
	///////////////////////////////////////////////////////////////////////////////////
	deque<OligoMatch> o_match;

	find_oligo_match(o_match, f_match, FORWARD, Seq_strand_plus, m_keys, m_db, m_seq);
	find_oligo_match(o_match, r_match, REVERSE, Seq_strand_minus, m_keys, m_db, m_seq);
		
	// Sort by sequence id, then sort by oligo location
	sort( o_match.begin(), o_match.end() );
		
	extract_amplicon_seq(amplicons, o_match, FORWARD, REVERSE, m_seq, m_amplicon_range,
		m_bounds_ptr);
	
	///////////////////////////////////////////////////////////////////////////////////
	// Find the set of sequences matched by {R(+), F(-)}
	///////////////////////////////////////////////////////////////////////////////////
	o_match.clear();

	find_oligo_match(o_match, f_match, FORWARD, Seq_strand_minus, m_keys, m_db, m_seq);
	find_oligo_match(o_match, r_match, REVERSE, Seq_strand_plus, m_keys, m_db, m_seq);

	// First sort by sequence id, then sort by oligo location
	sort( o_match.begin(), o_match.end() );
		
	extract_amplicon_seq(amplicons, o_match, REVERSE, FORWARD, m_seq, m_amplicon_range,
		m_bounds_ptr);
	
	// Make the amplicon sequences unique (since identical amplicons can be produced from different
	// target sequences).
	sort( amplicons.begin(), amplicons.end() );
	amplicons.erase( unique( amplicons.begin(), amplicons.end() ), amplicons.end() );
	
	// Convert each string into a Sequence
	return deque<Sequence>( amplicons.begin(), amplicons.end() );
}

bool PCR::multiplex_compatible(NucCruc &m_melt, 
	const Options &m_opt, const PCR &m_assay) const
{
	// To be conservative, we don't correct the oligo concentration for the oligo degeneracy
	const float effective_strand = m_opt.primer_strand;
	
	m_melt.strand(effective_strand);
	
	for(Oligo *query_oligo = assay_oligos();*query_oligo != LAST_OLIGO;++query_oligo){
	
		Word query_iter = oligo(*query_oligo).begin();

		do{ // Enumerate possibly degenerate oligos

			m_melt.set_query( query_iter.str() );

			for(Oligo *subject_oligo = m_assay.assay_oligos();
				*subject_oligo != LAST_OLIGO;++subject_oligo){
		
				// Test the query oligo against the subject oligo
				Word subject_iter = m_assay.oligo(*subject_oligo).begin();
				
				do{ // Enumerate possibly degenerate oligos

					m_melt.set_target( subject_iter.str() );

					if(m_melt.approximate_tm_heterodimer() >= m_opt.max_dimer){
						return false;
					}					

				} while( m_assay.oligo(*subject_oligo).next(subject_iter) );
			}

		} while( oligo(*query_oligo).next(query_iter) );
	}
	
	return true;
}
