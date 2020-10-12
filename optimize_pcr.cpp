#include <algorithm>
#include <math.h>
#include "assay.h"
#include "seq_overlap.h"

using namespace std;

pair<Word, Score> PCR::increase_degeneracy(const Oligo &m_oligo,
	const vector<Word> &m_target_keys, const float &m_target_threshold,
	const vector<Word> &m_background_keys, const float &m_background_threshold,
	const vector<Word> &m_multiplex_background_keys,
	const Score &m_score_threshold, NucCruc &m_melt, 
	const deque<PCR> &m_pool, const Options &m_opt)
{
	pair<Word, Score> ret;
	Score trial_score;
	
	if(oligo(m_oligo).degeneracy() >= m_opt.degen){
		return ret;
	}
	
	bool target_modified = false;
	bool background_modified = false;

	Word trial_oligo = oligo(m_oligo);

	// For multiplex assay design: the maximum overlap of all assay oligos
	// *except* m_oligo
	float partial_overlap = 0.0f;
	
	if(m_opt.use_multiplex){
		
		for(deque<PCR>::const_iterator i = m_pool.begin();i != m_pool.end();++i){
		
			// Since forward and reverse primers are interchangable, we need to 
			// compute the similarity to both primers. partial_overlap will store the
			// overlap of the oligo that is *not* m_oligo
			if(m_oligo == FORWARD){
			
				partial_overlap = max( partial_overlap, oligo(REVERSE).max_overlap( i->oligo(FORWARD) ) );
				partial_overlap = max( partial_overlap, oligo(REVERSE).max_overlap( i->oligo(REVERSE) ) );
			}
			else{ // m_oligo == REVERSE
				
				partial_overlap = max( partial_overlap, oligo(FORWARD).max_overlap( i->oligo(FORWARD) ) );
				partial_overlap = max( partial_overlap, oligo(FORWARD).max_overlap( i->oligo(REVERSE) ) );
			}
		}
		
		if(partial_overlap == 1.0f){
			partial_overlap = MULTIPLEX_OLIGO_REUSE_BONUS;
		}
	}
	
	const int first = trial_oligo.start();
	const int last = trial_oligo.stop();

	for(int i = first;i <= last;++i){

		for(unsigned char b = Base::A;b <= Base::T;b <<= 1){

			// Is the current base already degenerate and inclusive of b?
			if(trial_oligo.get(i) & b){
				continue;
			}
			
			// Increase the degeneracy by adding the b bit
			trial_oligo.mask(b, i);
			
			if( (trial_oligo.degeneracy() > m_opt.degen) || 
				!is_valid(m_oligo, trial_oligo, m_melt, m_opt, false /*no dimer check*/) ){
			
				// Restore the original base
				trial_oligo.unmask(b, i);

				continue;
			}

			switch(m_oligo){
				case FORWARD:
					update_identity(target_f_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
					break;
				case REVERSE:
					update_identity(target_r_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
					break;
				default:
					throw __FILE__ ":PCR::increase_degeneracy: Unknown oligo";
			};

			target_modified = true;

			trial_score.target_coverage = compute_target_coverage(m_opt.target_threshold);

			const float coverage_bound = trial_score.target_coverage + 
				m_score_threshold.background_coverage - 
				m_score_threshold.target_coverage;
			
			// If we are designing multiplex assays, we will keep testing an
			// assay with a coverage bound of zero (since it may yeild an improved
			// oligo overlap score).
			if( ( m_opt.use_multiplex && (coverage_bound < 0.0f) ) || 
			    ( !m_opt.use_multiplex && (coverage_bound <= 0.0f) ) ){

					// Restore the original base
					trial_oligo.unmask(b, i);

					continue; // Don't bother evaluating the background
			}

			switch(m_oligo){
				case FORWARD:
					update_identity(background_f_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
					update_identity(multiplex_background_f_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
					break;
				case REVERSE:
					update_identity(background_r_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
					update_identity(multiplex_background_r_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
					break;
				default:
					throw __FILE__ ":PCR::increase_degeneracy: Unknown oligo";
			};

			background_modified = true;

			trial_score.background_coverage = compute_background_coverage(m_opt.background_threshold);

			if(m_opt.use_multiplex){
				
				trial_score.background_coverage += compute_multiplex_background_coverage(m_opt.background_threshold);
				
				for(deque<PCR>::const_iterator m = m_pool.begin();m != m_pool.end();++m){
		
					// Since forward and reverse primers are interchangable, we need to 
					// compute the similarity to both primers
					trial_score.oligo_overlap = max( trial_score.oligo_overlap, 
						trial_oligo.max_overlap( m->oligo(FORWARD) ) );
					trial_score.oligo_overlap = max( trial_score.oligo_overlap, 
						trial_oligo.max_overlap( m->oligo(REVERSE) ) );
				}
				
				// Add the overlap for the oligos that we are not currenlty modifying
				trial_score.oligo_overlap = 
					( (trial_score.oligo_overlap == 1.0) ? MULTIPLEX_OLIGO_REUSE_BONUS : trial_score.oligo_overlap ) +
					partial_overlap;
			}
			
			if(trial_score > ret.second){

				ret.second = trial_score;
				ret.first = trial_oligo;
			}

			// Restore the original base
			trial_oligo.unmask(b, i);
		}			
	}

	if(target_modified){

		// Restore the indentity table
		switch(m_oligo){
			case FORWARD:
				update_identity(target_f_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(target_r_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::increase_degeneracy: Unknown oligo";
		};
	}

	if(background_modified){

		// Restore the indentity table
		switch(m_oligo){
			case FORWARD:
				update_identity(background_f_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_f_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(background_r_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_r_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::increase_degeneracy: Unknown oligo";
		};
	}
	
	return ret;
}

pair<Word, Score> PCR::decrease_degeneracy(const Oligo &m_oligo,
	const vector<Word> &m_target_keys, const float &m_target_threshold,
	const vector<Word> &m_background_keys, const float &m_background_threshold,
	const vector<Word> &m_multiplex_background_keys,
	const Score &m_score_threshold, NucCruc &m_melt, 
	const deque<PCR> &m_pool, const Options &m_opt)
{
	pair<Word, Score> ret;
	Word trial_oligo = oligo(m_oligo);
	Score trial_score;
	
	// For multiplex assay design: the maximum overlap of all assay oligos
	// *except* m_oligo
	float partial_overlap = 0.0f;
	
	if(m_opt.use_multiplex){
		
		for(deque<PCR>::const_iterator i = m_pool.begin();i != m_pool.end();++i){
		
			// Since forward and reverse primers are interchangable, we need to 
			// compute the similarity to both primers. partial_overlap will store the
			// overlap of the oligo that is *not* m_oligo
			if(m_oligo == FORWARD){
			
				partial_overlap = max( partial_overlap, oligo(REVERSE).max_overlap( i->oligo(FORWARD) ) );
				partial_overlap = max( partial_overlap, oligo(REVERSE).max_overlap( i->oligo(REVERSE) ) );
			}
			else{ // m_oligo == REVERSE
				
				partial_overlap = max( partial_overlap, oligo(FORWARD).max_overlap( i->oligo(FORWARD) ) );
				partial_overlap = max( partial_overlap, oligo(FORWARD).max_overlap( i->oligo(REVERSE) ) );
			}
		}
		
		if(partial_overlap == 1.0f){
			partial_overlap = MULTIPLEX_OLIGO_REUSE_BONUS;
		}
	}
	
	const int last = trial_oligo.stop();
	bool target_modified = false;
	bool background_modified = false;
		
	for(int i = trial_oligo.start();i <= last;++i){
		
		const unsigned char curr = trial_oligo.get(i);
		
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
			trial_oligo.unmask(b, i);
			
			if( !is_valid(m_oligo, trial_oligo, m_melt, m_opt, false /*no dimer check*/) ){
			
				// Restore the original base
				trial_oligo.mask(b, i);

				continue;
			}
			
			switch(m_oligo){
				case FORWARD:
					update_identity(target_f_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
					break;
				case REVERSE:
					update_identity(target_r_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
					break;
				default:
					throw __FILE__ ":PCR::decrease_degeneracy: Unknown oligo";
			};
			
			target_modified = true;
			
			trial_score.target_coverage = compute_target_coverage(m_opt.target_threshold);
			
			const float coverage_bound = trial_score.target_coverage + 
				m_score_threshold.background_coverage - 
				m_score_threshold.target_coverage;
			
			// If we are designing multiplex assays, we will keep testing an
			// assay with a coverage bound of zero (since it may yeild an improved
			// oligo overlap score).
			if( ( m_opt.use_multiplex && (coverage_bound < 0.0f) ) || 
			    ( !m_opt.use_multiplex && (coverage_bound <= 0.0f) ) ){
					
					// Restore the original base
					trial_oligo.mask(b, i);
					
					continue; // Don't bother evaluating the background
			}
					
			switch(m_oligo){
				case FORWARD:
					update_identity(background_f_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
					update_identity(multiplex_background_f_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
					break;
				case REVERSE:
					update_identity(background_r_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
					update_identity(multiplex_background_r_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
					break;
				default:
					throw __FILE__ ":decrease_degeneracy: Unknown oligo";
			};
			
			background_modified = true;
			
			trial_score.background_coverage = compute_background_coverage(m_opt.background_threshold);

			if(m_opt.use_multiplex){
				
				trial_score.background_coverage += compute_multiplex_background_coverage(m_opt.background_threshold);
				
				trial_score.oligo_overlap = 0.0f;
				
				for(deque<PCR>::const_iterator m = m_pool.begin();m != m_pool.end();++m){
		
					// Since forward and reverse primers are interchangable, we need to 
					// compute the similarity to both primers
					trial_score.oligo_overlap = max( trial_score.oligo_overlap, 
						trial_oligo.max_overlap( m->oligo(FORWARD) ) );
					trial_score.oligo_overlap = max( trial_score.oligo_overlap, 
						trial_oligo.max_overlap( m->oligo(REVERSE) ) );
				}
				
				// Add the overlap for the oligos that we are not currenlty modifying
				trial_score.oligo_overlap = 
					( (trial_score.oligo_overlap == 1.0f) ? MULTIPLEX_OLIGO_REUSE_BONUS : trial_score.oligo_overlap ) +
					partial_overlap;
			}
			
			if(trial_score > ret.second){

				ret.second = trial_score;
				ret.first = trial_oligo;
			}
			
			// Restore the original base
			trial_oligo.mask(b, i);
		}			
	}

	if(target_modified){
		
		// Restore the indentity table
		switch(m_oligo){
			case FORWARD:
				update_identity(target_f_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(target_r_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":increase_degeneracy: Unknown oligo";
		};
	}

	if(background_modified){

		// Restore the indentity table
		switch(m_oligo){
			case FORWARD:
				update_identity(background_f_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_f_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(background_r_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_r_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":increase_degeneracy: Unknown oligo";
		};
	}
	
	return ret;
}

pair<Word, Score> PCR::trim5(const Oligo &m_oligo,
	const vector<Word> &m_target_keys, const float &m_target_threshold,
	const vector<Word> &m_background_keys, const float &m_background_threshold,
	const vector<Word> &m_multiplex_background_keys,
	const Score &m_score_threshold, NucCruc &m_melt, 
	const deque<PCR> &m_pool, const Options &m_opt)
{
	const int len = oligo(m_oligo).size();
	pair<Word, Score> ret;
	Score trial_score;
		
	// Already at the minimum length, don't shrink further
	if(len == m_opt.primer_range.first){
		return ret;
	}

	Word trial_oligo = oligo(m_oligo);
	
	// Remove the 5' base
	trial_oligo.shrink_front(); // <-- set the first valid base to EOS
	
	if( !is_valid(m_oligo, trial_oligo, m_melt, m_opt, false /*no dimer check*/) ){
		return ret;
	}
			
	switch(m_oligo){
		case FORWARD:
			update_identity(target_f_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
			break;
		case REVERSE:
			update_identity(target_r_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
			break;
		default:
			throw __FILE__ ":PCR::trim5: Unknown oligo";
	};

	trial_score.target_coverage = compute_target_coverage(m_opt.target_threshold);
			
	const float coverage_bound = trial_score.target_coverage + 
		m_score_threshold.background_coverage - 
		m_score_threshold.target_coverage;

	// If we are designing multiplex assays, we will keep testing an
	// assay with a coverage bound of zero (since it may yeild an improved
	// oligo overlap score).
	if( ( m_opt.use_multiplex && (coverage_bound < 0.0f) ) || 
	    ( !m_opt.use_multiplex && (coverage_bound <= 0.0f) ) ){
			
			// Restore the target indentity table
			switch(m_oligo){
				case FORWARD:
					update_identity(target_f_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
					break;
				case REVERSE:
					update_identity(target_r_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
					break;
				default:
					throw __FILE__ ":PCR::trim5: Unknown oligo";
			};
			
			return ret; // Don't bother evaluating the background
	}
					
	switch(m_oligo){
		case FORWARD:
			update_identity(background_f_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
			update_identity(multiplex_background_f_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
			break;
		case REVERSE:
			update_identity(background_r_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
			update_identity(multiplex_background_r_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
			break;
		default:
			throw __FILE__ ":PCR::trim5: Unknown oligo";
	};
	
	trial_score.background_coverage = compute_background_coverage(m_opt.background_threshold);
	
	if(m_opt.use_multiplex){
		
		trial_score.background_coverage += compute_multiplex_background_coverage(m_opt.background_threshold);
		
		float best_f = 0.0;
		float best_r = 0.0;
		
		for(deque<PCR>::const_iterator i = m_pool.begin();i != m_pool.end();++i){
		
			// Since forward and reverse primers are interchangable, we need to 
			// compute the similarity to both primers
			if(m_oligo == FORWARD){
			
				best_f = max( best_f, trial_oligo.max_overlap( i->oligo(FORWARD) ) );
				best_f = max( best_f, trial_oligo.max_overlap( i->oligo(REVERSE) ) );
			}
			else{
				best_f = max( best_f, oligo(FORWARD).max_overlap( i->oligo(FORWARD) ) );
				best_f = max( best_f, oligo(FORWARD).max_overlap( i->oligo(REVERSE) ) );
			}
			
			if(m_oligo == REVERSE){
			
				best_r = max( best_r, trial_oligo.max_overlap( i->oligo(FORWARD) ) );
				best_r = max( best_r, trial_oligo.max_overlap( i->oligo(REVERSE) ) );
			}
			else{
				best_r = max( best_r, oligo(REVERSE).max_overlap( i->oligo(FORWARD) ) );
				best_r = max( best_r, oligo(REVERSE).max_overlap( i->oligo(REVERSE) ) );
			}
		}
		
		trial_score.oligo_overlap = 
					( (best_f == 1.0f) ? MULTIPLEX_OLIGO_REUSE_BONUS : best_f ) +
					( (best_r == 1.0f) ? MULTIPLEX_OLIGO_REUSE_BONUS : best_r );
		
	}
	
	if(trial_score > ret.second){

		ret.second = trial_score;
		ret.first = trial_oligo;
	}
	
	// Restore the target and background indentity tables
	switch(m_oligo){
		case FORWARD:
			update_identity(target_f_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
			update_identity(background_f_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
			update_identity(multiplex_background_f_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
			break;
		case REVERSE:
			update_identity(target_r_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
			update_identity(background_r_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
			update_identity(multiplex_background_r_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
			break;
		default:
			throw __FILE__ ":PCR::trim5: Unknown oligo";
	};
	
	return ret;
}

pair<Word, Score> PCR::trim3(const Oligo &m_oligo,
	const vector<Word> &m_target_keys, const float &m_target_threshold,
	const vector<Word> &m_background_keys, const float &m_background_threshold,
	const vector<Word> &m_multiplex_background_keys,
	const Score &m_score_threshold, NucCruc &m_melt, 
	const deque<PCR> &m_pool, const Options &m_opt)
{
	const int len = oligo(m_oligo).size();
	pair<Word, Score> ret;
	Score trial_score;
	
	// Already at the minimum length, don't shrink further
	if(len == m_opt.primer_range.first){
		return ret;
	}

	Word trial_oligo = oligo(m_oligo);
	
	// Remove the 3' base
	trial_oligo.shrink_back(); // <-- set the last valid base to EOS
	
	if( !is_valid(m_oligo, trial_oligo, m_melt, m_opt, false /*no dimer check*/) ){
		return ret;
	}
	
	switch(m_oligo){
		case FORWARD:
			update_identity(target_f_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
			break;
		case REVERSE:
			update_identity(target_r_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
			break;
		default:
			throw __FILE__ ":PCR::trim3: Unknown oligo";
	};

	trial_score.target_coverage = compute_target_coverage(m_opt.target_threshold);
	
	const float coverage_bound = trial_score.target_coverage + 
		m_score_threshold.background_coverage - 
		m_score_threshold.target_coverage;

	// If we are designing multiplex assays, we will keep testing an
	// assay with a coverage bound of zero (since it may yeild an improved
	// oligo overlap score).
	if( ( m_opt.use_multiplex && (coverage_bound < 0.0f) ) || 
	    ( !m_opt.use_multiplex && (coverage_bound <= 0.0f) ) ){
			
			// Restore the target indentity table
			switch(m_oligo){
				case FORWARD:
					update_identity(target_f_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
					break;
				case REVERSE:
					update_identity(target_r_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
					break;
				default:
					throw __FILE__ ":PCR::trim5: Unknown oligo";
			};
			
			return ret; // Don't bother evaluating the background
	}
	
	switch(m_oligo){
		case FORWARD:
			update_identity(background_f_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
			update_identity(multiplex_background_f_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
			break;
		case REVERSE:
			update_identity(background_r_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
			update_identity(multiplex_background_r_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
			break;
		default:
			throw __FILE__ ":PCR::trim3: Unknown oligo";
	};
	
	trial_score.background_coverage = compute_background_coverage(m_opt.background_threshold);
	
	if(m_opt.use_multiplex){
		
		trial_score.background_coverage += compute_multiplex_background_coverage(m_opt.background_threshold);
		
		float best_f = 0.0;
		float best_r = 0.0;
		
		for(deque<PCR>::const_iterator i = m_pool.begin();i != m_pool.end();++i){
		
			// Since forward and reverse primers are interchangable, we need to 
			// compute the similarity to both primers
			if(m_oligo == FORWARD){
			
				best_f = max( best_f, trial_oligo.max_overlap( i->oligo(FORWARD) ) );
				best_f = max( best_f, trial_oligo.max_overlap( i->oligo(REVERSE) ) );
			}
			else{
				best_f = max( best_f, oligo(FORWARD).max_overlap( i->oligo(FORWARD) ) );
				best_f = max( best_f, oligo(FORWARD).max_overlap( i->oligo(REVERSE) ) );
			}
			
			if(m_oligo == REVERSE){
			
				best_r = max( best_r, trial_oligo.max_overlap( i->oligo(FORWARD) ) );
				best_r = max( best_r, trial_oligo.max_overlap( i->oligo(REVERSE) ) );
			}
			else{
				best_r = max( best_r, oligo(REVERSE).max_overlap( i->oligo(FORWARD) ) );
				best_r = max( best_r, oligo(REVERSE).max_overlap( i->oligo(REVERSE) ) );
			}
		}
		
		trial_score.oligo_overlap = 
					( (best_f == 1.0f) ? MULTIPLEX_OLIGO_REUSE_BONUS : best_f ) +
					( (best_r == 1.0f) ? MULTIPLEX_OLIGO_REUSE_BONUS : best_r );
	}
	
	if(trial_score > ret.second){

		ret.second = trial_score;
		ret.first = trial_oligo;
	}

	// Restore the indentity tables
	switch(m_oligo){
		case FORWARD:
			update_identity(target_f_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
			update_identity(background_f_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
			update_identity(multiplex_background_f_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
			break;
		case REVERSE:
			update_identity(target_r_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
			update_identity(background_r_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
			update_identity(multiplex_background_r_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
			break;
		default:
			throw __FILE__ ":PCR::trim3: Unknown oligo";
	};
		
	return ret;
}

pair<Word, Score> PCR::grow5(const Oligo &m_oligo,
	const vector<Word> &m_target_keys, const float &m_target_threshold,
	const vector<Word> &m_background_keys, const float &m_background_threshold,
	const vector<Word> &m_multiplex_background_keys,
	const Score &m_score_threshold, NucCruc &m_melt, 
	const deque<PCR> &m_pool, const Options &m_opt)
{
	const int len = oligo(m_oligo).size();
	pair<Word, Score> ret;
	Score trial_score;
	
	// Already at the maximum length, don't grow further
	if(len == m_opt.primer_range.second){
		return ret;
	}
	
	bool target_modified = false;
	bool background_modified = false;
	
	// For multiplex assay design: the maximum overlap of all assay oligos
	// *except* m_oligo
	float partial_overlap = 0.0f;
	
	if(m_opt.use_multiplex){
	
		for(deque<PCR>::const_iterator i = m_pool.begin();i != m_pool.end();++i){

			// Since forward and reverse primers are interchangable, we need to 
			// compute the similarity to both primers. partial_overlap will store the
			// overlap of the oligo that is *not* m_oligo
			if(m_oligo == FORWARD){

				partial_overlap = max( partial_overlap, oligo(REVERSE).max_overlap( i->oligo(FORWARD) ) );
				partial_overlap = max( partial_overlap, oligo(REVERSE).max_overlap( i->oligo(REVERSE) ) );
			}
			else{ // m_oligo == REVERSE

				partial_overlap = max( partial_overlap, oligo(FORWARD).max_overlap( i->oligo(FORWARD) ) );
				partial_overlap = max( partial_overlap, oligo(FORWARD).max_overlap( i->oligo(REVERSE) ) );
			}
		}
		
		if(partial_overlap == 1.0f){
			partial_overlap = MULTIPLEX_OLIGO_REUSE_BONUS;
		}
	}
	
	for(unsigned char b = Base::A;b <= Base::T;b <<= 1){

		Word trial_oligo = oligo(m_oligo);
		
		trial_oligo.grow_front(b); // If there's room, add b to the front of the oligo
		
		if( !is_valid(m_oligo, trial_oligo, m_melt, m_opt, false /*no dimer check*/) ){
			continue;
		}
	
		switch(m_oligo){
			case FORWARD:
				update_identity(target_f_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(target_r_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::grow5: Unknown oligo";
		};
		
		target_modified = true;
		
		trial_score.target_coverage = compute_target_coverage(m_opt.target_threshold);
		
		const float coverage_bound = trial_score.target_coverage + 
			m_score_threshold.background_coverage - 
			m_score_threshold.target_coverage;

		// If we are designing multiplex assays, we will keep testing an
		// assay with a coverage bound of zero (since it may yeild an improved
		// oligo overlap score).
		if( ( m_opt.use_multiplex && (coverage_bound < 0.0f) ) || 
		    ( !m_opt.use_multiplex && (coverage_bound <= 0.0f) ) ){
				continue; // Don't bother evaluating the background
		}

		switch(m_oligo){
			case FORWARD:
				update_identity(background_f_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_f_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(background_r_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_r_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::grow5: Unknown oligo";
		};
		
		background_modified = true;
		
		trial_score.background_coverage = compute_background_coverage(m_opt.background_threshold);

		if(m_opt.use_multiplex){
			
			trial_score.background_coverage += compute_multiplex_background_coverage(m_opt.background_threshold);
			
			trial_score.oligo_overlap = 0.0f;
			
			for(deque<PCR>::const_iterator m = m_pool.begin();m != m_pool.end();++m){

				// Since forward and reverse primers are interchangable, we need to 
				// compute the similarity to both primers
				trial_score.oligo_overlap = max( trial_score.oligo_overlap, 
					trial_oligo.max_overlap( m->oligo(FORWARD) ) );
				trial_score.oligo_overlap = max( trial_score.oligo_overlap, 
					trial_oligo.max_overlap( m->oligo(REVERSE) ) );
			}

			// Add the overlap for the oligos that we are not currenlty modifying
			trial_score.oligo_overlap = 
					( (trial_score.oligo_overlap == 1.0f) ? MULTIPLEX_OLIGO_REUSE_BONUS : trial_score.oligo_overlap ) +
					partial_overlap;
		}
		
		if(trial_score > ret.second){

			ret.second = trial_score;
			ret.first = trial_oligo;
		}
	}

	if(target_modified){
	
		// Restore the indentity table
		switch(m_oligo){
			case FORWARD:
				update_identity(target_f_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(target_r_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::grow5: Unknown oligo";
		};
	}
	
	if(background_modified){
	
		// Restore the indentity table
		switch(m_oligo){
			case FORWARD:
				update_identity(background_f_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_f_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(background_r_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_r_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::grow5: Unknown oligo";
		};
	}
	
	return ret;
}

pair<Word, Score> PCR::grow3(const Oligo &m_oligo,
	const vector<Word> &m_target_keys, const float &m_target_threshold,
	const vector<Word> &m_background_keys, const float &m_background_threshold,
	const vector<Word> &m_multiplex_background_keys,
	const Score &m_score_threshold, NucCruc &m_melt, 
	const deque<PCR> &m_pool, const Options &m_opt)
{
	pair<Word, Score> ret;
	Score trial_score;
	
	const int len = oligo(m_oligo).size();
	
	// Already at the maximum length, don't grow further
	if(len == m_opt.primer_range.second){
		return ret;
	}
	
	// For multiplex assay design: the maximum overlap of all assay oligos
	// *except* m_oligo
	float partial_overlap = 0.0f;
	
	if(m_opt.use_multiplex){
	
		for(deque<PCR>::const_iterator i = m_pool.begin();i != m_pool.end();++i){

			// Since forward and reverse primers are interchangable, we need to 
			// compute the similarity to both primers. partial_overlap will store the
			// overlap of the oligo that is *not* m_oligo
			if(m_oligo == FORWARD){

				partial_overlap = max( partial_overlap, oligo(REVERSE).max_overlap( i->oligo(FORWARD) ) );
				partial_overlap = max( partial_overlap, oligo(REVERSE).max_overlap( i->oligo(REVERSE) ) );
			}
			else{ // m_oligo == REVERSE

				partial_overlap = max( partial_overlap, oligo(FORWARD).max_overlap( i->oligo(FORWARD) ) );
				partial_overlap = max( partial_overlap, oligo(FORWARD).max_overlap( i->oligo(REVERSE) ) );
			}
		}
		
		if(partial_overlap == 1.0f){
			partial_overlap = MULTIPLEX_OLIGO_REUSE_BONUS;
		}
	}
	
	bool target_modified = false;
	bool background_modified = false;
	
	for(unsigned char b = Base::A;b <= Base::T;b <<= 1){
		
		Word trial_oligo = oligo(m_oligo);
		
		trial_oligo.grow_back(b); // If there's room, add b to the back of the oligo
		
		if( !is_valid(m_oligo, trial_oligo, m_melt, m_opt, false /*no dimer check*/) ){
			continue;
		}
		
		switch(m_oligo){
			case FORWARD:
				update_identity(target_f_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(target_r_identity, trial_oligo, m_target_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::grow3: Unknown oligo";
		};
		
		target_modified = true;
		
		trial_score.target_coverage = compute_target_coverage(m_opt.target_threshold);
		
		const float coverage_bound = trial_score.target_coverage + 
			m_score_threshold.background_coverage - 
			m_score_threshold.target_coverage;

		// If we are designing multiplex assays, we will keep testing an
		// assay with a coverage bound of zero (since it may yeild an improved
		// oligo overlap score).
		if( ( m_opt.use_multiplex && (coverage_bound < 0.0f) ) || 
		    ( !m_opt.use_multiplex && (coverage_bound <= 0.0f) ) ){
				continue; // Don't bother evaluating the background
		}

		switch(m_oligo){
			case FORWARD:
				update_identity(background_f_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_f_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(background_r_identity, trial_oligo, m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_r_identity, trial_oligo, m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::grow3: Unknown oligo";
		};
		
		background_modified = true;
		
		trial_score.background_coverage = compute_background_coverage(m_opt.background_threshold);

		if(m_opt.use_multiplex){
			
			trial_score.background_coverage += compute_multiplex_background_coverage(m_opt.background_threshold);
			
			trial_score.oligo_overlap = 0.0f;
			
			for(deque<PCR>::const_iterator m = m_pool.begin();m != m_pool.end();++m){

				// Since forward and reverse primers are interchangable, we need to 
				// compute the similarity to both primers
				trial_score.oligo_overlap = max( trial_score.oligo_overlap, 
					trial_oligo.max_overlap( m->oligo(FORWARD) ) );
				trial_score.oligo_overlap = max( trial_score.oligo_overlap, 
					trial_oligo.max_overlap( m->oligo(REVERSE) ) );
			}

			// Add the overlap for the oligos that we are not currenlty modifying
			trial_score.oligo_overlap = 
					( (trial_score.oligo_overlap == 1.0f) ? MULTIPLEX_OLIGO_REUSE_BONUS : trial_score.oligo_overlap ) +
					partial_overlap;
		}
		
		if(trial_score > ret.second){

			ret.second = trial_score;
			ret.first = trial_oligo;
		}
	}

	if(target_modified){
		
		// Restore the indentity table
		switch(m_oligo){
			case FORWARD:
				update_identity(target_f_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(target_r_identity, oligo(m_oligo), m_target_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::grow3: Unknown oligo";
		};
	}

	if(background_modified){

		// Restore the indentity table
		switch(m_oligo){
			case FORWARD:
				update_identity(background_f_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_f_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			case REVERSE:
				update_identity(background_r_identity, oligo(m_oligo), m_background_keys, m_opt.use_taq_mama);
				update_identity(multiplex_background_r_identity, oligo(m_oligo), m_multiplex_background_keys, m_opt.use_taq_mama);
				break;
			default:
				throw __FILE__ ":PCR::grow3: Unknown oligo";
		};
	}
	
	return ret;
}
