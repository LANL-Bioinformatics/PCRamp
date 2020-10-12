#include <algorithm>
#include <unordered_set>
#include <math.h>
#include "assay.h"
#include "seq_overlap.h"

using namespace std;

const char* oligo_name[LAST_OLIGO] = {
	"F",
	"R"
};

Score optimize(PCR &m_assay,
	const vector<Move> &m_moves,
	const vector<Word> &m_target_keys,
	const MULTIMAP<Word, WordMatch> &m_target_db,
	const deque<Sequence> &m_target_seq,
	
	const vector<Word> &m_background_keys,
	const MULTIMAP<Word, WordMatch> &m_background_db,
	const deque<Sequence> &m_background_seq,
	
	const vector<Word> &m_multiplex_background_keys,
	const MULTIMAP<Word, WordMatch> &m_multiplex_background_db,
	const deque<Sequence> &m_multiplex_background_seq,
	
	const deque<PCR> &m_pool,
	const Options &m_opt,
	ostream &m_out)
{	
	
	PCR best;
	Score best_score;
		
	PCR approx;
	Score approx_score;
	
	best.copy_oligos(m_assay);
	approx.copy_oligos(m_assay);
	
	unsigned int iteration = 0;
	
	// Make sure not to get stuck in a loop (endlessly growing and shrinking the same assays)
	unordered_set<string> previous;
	
	previous.insert( best.packed_string() );
	
	NucCruc melt;
	
	melt.fast_alignment(true); // No gaps allowed
	melt.salt(m_opt.salt);
	
	while(true){
		
		bool improved = false;
		
		++iteration;
		
		// Store the potential assay match candidates for each sequence
		approx.collect_target_candidates(
			m_target_keys, m_target_db, m_target_seq,
			m_opt);
		
		approx.collect_background_candidates(
			m_background_keys, m_background_db, m_background_seq,
			m_opt);
				
		approx.update_target_candidates(m_target_keys, m_opt.use_taq_mama);
		
		approx.update_background_candidates(m_background_keys, m_opt.use_taq_mama);
		
		approx_score.target_coverage = 
			approx.compute_target_coverage(m_opt.target_threshold);
		
		approx_score.background_coverage = 
			approx.compute_background_coverage(m_opt.background_threshold);
				
		if(m_opt.use_multiplex){
			
			// If we have existing amplicon sequences for PCR-based assays, we
			// treat them like backgrounds
			approx.collect_multiplex_background_candidates(
				m_multiplex_background_keys, m_multiplex_background_db, m_multiplex_background_seq,
				m_opt);
			
			approx.update_multiplex_background_candidates(m_multiplex_background_keys, m_opt.use_taq_mama);
			
			// Add the multiplex background coverage to the existing background background
			approx_score.background_coverage += 
				approx.compute_multiplex_background_coverage(m_opt.background_threshold);
			
			// When multiplexing assays, encourage oligo reuse by optimizing the
			// amount of similarity with previously designed oligos
			approx_score.oligo_overlap = 
				approx.compute_oligo_overlap(m_pool);
		}
			
		if( (iteration == 1) && (m_opt.output_filter > Options::VERBOSE) ){
			m_out << "\t\tinitial accuracy = " << approx_score.accuracy() << endl;
		}
		
		if(approx_score < best_score){
			
			// No improvement -- we can get here because the optimization algorithm
			// may change some assay parameters (like amplicon length) so that
			// they are no longer valid
			break;
		}
		
		best_score = approx_score;
		best.copy_oligos(approx);
		
		Word local_seq;
		Oligo local_oligo = LAST_OLIGO;
		Move local_move = LastMove;
		Score local_score = approx_score;
		
		//////////////////////////////////////////////////////////////////////
		// Search each assay oligo
		//////////////////////////////////////////////////////////////////////
		for(Oligo *o_ptr = approx.assay_oligos();*o_ptr != LAST_OLIGO;++o_ptr){
		
			for(vector<Move>::const_iterator m = m_moves.begin();m != m_moves.end();++m){

				pair<Word, Score> tmp = optimization_move(*m, *o_ptr, approx, 
					m_target_keys, m_opt.target_threshold,
					m_background_keys, m_opt.background_threshold, 
					m_multiplex_background_keys,
					local_score, melt, m_pool, m_opt);
				
				// Select for the highest score, or, if scores are equal, the lowest degeneracy
				if( (tmp.second > local_score) || 
					( (tmp.second == local_score) && 
					  (tmp.first.degeneracy() < local_seq.degeneracy() ) ) ){
			
					local_score = tmp.second;
					local_seq = tmp.first;
					local_oligo = *o_ptr;
					local_move = *m;
					improved = true;
				}
			}
		}

		if(!improved){
			break;
		}
		
		approx_score = local_score;
		
		// Center the modified oligo in the word
		local_seq.center();
		
		approx.oligo(local_oligo, local_seq);
		
		if(m_opt.output_filter > Options::VERBOSE){
		
			m_out << "\t\tapprox accuracy[" << iteration << "] = " << approx_score.accuracy()
				<< " (" << approx_score.target_coverage << ", " << approx_score.background_coverage << ")" ;

			if(m_opt.use_multiplex){
				m_out << ":" << approx_score.oligo_overlap;
			}
			
			switch(local_move){
				case IncreaseDegeneracy:
					m_out << ":" << oligo_name[local_oligo] << " +Degen: ";
					break;
				case DecreaseDegeneracy:
					m_out << ":" << oligo_name[local_oligo] << " -Degen: ";
					break;
				case Trim5:
					m_out << ":" << oligo_name[local_oligo] << " -5': ";
					break;
				case Trim3: 
					m_out << ":" << oligo_name[local_oligo] << " -3': ";
					break;
				case Grow5: 
					m_out << ":" << oligo_name[local_oligo] << " +5': ";
					break;
				case Grow3: 
					m_out << ":" << oligo_name[local_oligo] << " +3': ";
					break;
				default:
					throw __FILE__ ":optimize: Unknown optimization move!";
			};

			approx.write(m_out);

			m_out << endl;
		}
		
		// Have we already tried this assay?
		const string p_str = approx.packed_string();
		
		if( previous.find(p_str) != previous.end() ){
			break;
		}
		
		previous.insert(p_str);
	}
	
	m_assay.copy_oligos(best);
			
	return best_score;
}

void update_identity(MAP<unsigned int, float> &m_ident, const Word &m_w, const vector<Word> &m_keys, bool m_use_taq_mama)
{
	typedef MAP<unsigned int, float>::iterator I;
	
	if( m_ident.empty() ){
		return;
	}
	
	const unsigned int len = m_w.size();
	
	assert(len > 0);
	
	const float norm = 1.0/len;
	
	for(I i = m_ident.begin();i != m_ident.end();++i){
		i->second = (m_w & m_keys[i->first])*norm;
	}
		
	if(!m_use_taq_mama){
		return;
	}
	
	// Apply the TaqMAMA correction
	const int last = m_w.stop();
	const int penultimate = last - 1;

	const pair<unsigned char, unsigned char> primer( m_w.get(penultimate), m_w.get(last) );

	// The Taq-MAMA correction is not defined for degenerate bases -- fall through to the
	// standard identity calculation if either of the last two primer bases are degenerate
	if( !is_degen(primer.first) && !is_degen(primer.second) ){

		for(I i = m_ident.begin();i != m_ident.end();++i){

			const Word &key_ref = m_keys[i->first];

			const pair<unsigned char, unsigned char> template_seq( key_ref.get(penultimate), key_ref.get(last) );

			if( !is_degen(template_seq.first) && !is_degen(template_seq.second) ){


				// When computing total number of mismatches, don't include
				// mismatches that occur at the last two 3' positions -- these
				// will be accounted for in the Taq-MAMA correction
				//i->second = ( (m_w & key_ref) + // total number of matches
				//	(primer.first != template_seq.first) + // add back a penultimate mismatch
				//	(primer.second != template_seq.second) )* // add back an ultimate mismatch
				//	norm*taq_mama_correction(primer, template_seq);
				i->second *= taq_mama_correction(primer, template_seq);
			}
		}
	}
}

void find_oligo_match(deque<OligoMatch> &m_match, const deque<unsigned int> &m_word_matches, const Oligo &m_oligo, 
	const Strand &m_strand, const vector<Word> &m_keys,
	const MULTIMAP<Word, WordMatch> &m_db,
	const deque<Sequence> &m_seq)
{
	typedef MULTIMAP<Word, WordMatch>::const_iterator I;
	
	for(deque<unsigned int>::const_iterator i = m_word_matches.begin();i != m_word_matches.end();++i){
		
		const pair<I, I> range = m_db.equal_range(m_keys[*i]);
		
		for(I j = range.first;j != range.second;++j){
			
			// Skip oligos that don't match the specified strand
			if( !(j->second.s & m_strand) ){
				continue;
			}
			
			// Skip sequences that are not active
			if( !m_seq[j->second.index].active() ){
				continue;
			}
			
			m_match.push_back( OligoMatch(m_oligo, *i, j->second) );
		}
	}
}
	
void match_words(deque<unsigned int> &m_match, const Word &m_seq, const vector<Word> &m_keys, const float &m_threshold)
{
	const unsigned int scaled_threshold = m_seq.size()*m_threshold;
		
	for(vector<Word>::const_iterator k = m_keys.begin();k != m_keys.end();++k){

		if( (m_seq & *k) >= scaled_threshold ){
			m_match.push_back( k - m_keys.begin() );
		}
	}
}

pair<Word, Score> optimization_move(const Move &m_move,
	const Oligo &m_oligo, 
	PCR &m_assay, 
	const vector<Word> &m_target_keys, const float &m_target_threshold,
	const vector<Word> &m_background_keys, const float &m_background_threshold,
	const vector<Word> &m_multiplex_background_keys,
	const Score &m_score_threshold, NucCruc &m_melt, 
	const deque<PCR> &m_pool, const Options &m_opt)
{	
	switch(m_move){
		case IncreaseDegeneracy:
			return m_assay.increase_degeneracy(m_oligo, 
				m_target_keys, m_target_threshold,
				m_background_keys, m_background_threshold,
				m_multiplex_background_keys,
				m_score_threshold, m_melt, m_pool, m_opt);
		case DecreaseDegeneracy:
			return m_assay.decrease_degeneracy(m_oligo, 
				m_target_keys, m_target_threshold,
				m_background_keys, m_background_threshold,
				m_multiplex_background_keys,
				m_score_threshold, m_melt, m_pool, m_opt);
		case Trim5:
			return m_assay.trim5(m_oligo, 
				m_target_keys, m_target_threshold,
				m_background_keys, m_background_threshold,
				m_multiplex_background_keys,
				m_score_threshold, m_melt, m_pool, m_opt);
		case Trim3:
			return m_assay.trim3(m_oligo, 
				m_target_keys, m_target_threshold,
				m_background_keys, m_background_threshold,
				m_multiplex_background_keys,
				m_score_threshold, m_melt, m_pool, m_opt);
		case Grow5:
			return m_assay.grow5(m_oligo, 
				m_target_keys, m_target_threshold,
				m_background_keys, m_background_threshold,
				m_multiplex_background_keys,
				m_score_threshold, m_melt, m_pool, m_opt);
		case Grow3:
			return m_assay.grow3(m_oligo, 
				m_target_keys, m_target_threshold,
				m_background_keys, m_background_threshold,
				m_multiplex_background_keys,
				m_score_threshold, m_melt, m_pool, m_opt);
		default:
			throw __FILE__ ":optimization_move: Unknown move";
	};
	
	return make_pair( Word(), Score() );
}

bool make_degenerate(PCR &m_assay, 
		const vector<Word> &m_target_keys, 
		const MULTIMAP<Word, WordMatch> &m_target_db, 
		const deque<Sequence> &m_target_seq,
		NucCruc &m_melt,  
		const Options &m_opt,
		ostream &m_out)
{
	PCR local;

	local.copy_oligos(m_assay);
	
	// Store the potential assay match candidates for each sequence
	local.collect_target_candidates(
		m_target_keys, m_target_db, m_target_seq,
		m_opt);

	// Our challenge is to maximize the number of targets covered (or at least maximize the degeneracy?)
	// while satisfying the thermodynamic constraints on each oligo (i.e. target Tm within allowed range,
	// no primer hairpin or primer dimer formation, etc.).
	local.update_target_candidates(m_target_keys, m_opt.use_taq_mama);

	local.sort_target_candidates();
	
	const bool ret = local.maximize_degeneracy(m_target_keys, m_melt, m_opt);
	
	m_assay.copy_oligos(local);
		
	if(m_opt.output_filter > Options::VERBOSE){

		#pragma omp critical
		if(ret){
			m_out << "Init top-down assay: ";
			m_assay.write(m_out);
			m_out << endl;
		}
		else{
			m_out << "Unable to initialize a top-down assay!" << endl;
		}
	}
	
	return ret;
}
