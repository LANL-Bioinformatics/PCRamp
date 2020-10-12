#include <math.h>
#include "assay.h"
#include "seq_overlap.h"

using namespace std;

void PCR::find_background_match(BitSet &m_match, const std::vector<Word> &m_keys, 
	const MULTIMAP<Word, WordMatch> &m_db, const deque<Sequence> &m_seq, 
	const Options &m_opt, std::ostream &m_out)
{
	const unsigned int num_seq = m_seq.size();
	const float perfect_match_score = PERFECT_MATCH_SCORE;
	
	if(m_match.size() != num_seq){
		throw __FILE__ ":PCR::find_background_match: |match| != num_seq";
	}
	
	collect_background_candidates(m_keys, m_db, m_seq, m_opt);
	
	float f_norm = perfect_match_score*oligo(FORWARD).size();
	float r_norm = perfect_match_score*oligo(REVERSE).size();
	
	if(f_norm > 0.0f){
		f_norm = 1.0f/f_norm;
	}

	if(r_norm > 0.0f){
		r_norm = 1.0f/r_norm;
	}
	
	pair<unsigned char, unsigned char> Fp_taq_mama;
	pair<unsigned char, unsigned char> Fm_taq_mama;
	pair<unsigned char, unsigned char> Rp_taq_mama;
	pair<unsigned char, unsigned char> Rm_taq_mama;
	
	if(m_opt.use_taq_mama){
		
		Fp_taq_mama = oligo(FORWARD).get_last_two();
		Fm_taq_mama = oligo(FORWARD).complement().get_last_two();
		Rp_taq_mama = oligo(REVERSE).get_last_two();
		Rm_taq_mama = oligo(REVERSE).complement().get_last_two();
	}
	
	SO::SeqOverlap align(SO::SeqOverlap::SmithWaterman, true /* is nucleic acid */);

	// Packing scheme:
	// slot 0: F   + f[i]
	// slot 1: (F) + f[i]
	// slot 2: R   + r[i]
	// slot 3: (R) + r[i]
	// slot 4: F   + f[i + 1]
	// slot 5: (F) + f[i + 1]
	// slot 6: R   + r[i + 1]
	// slot 7: (R) + r[i + 1]

	// Simple method #1. Search forward and reverse sequences in parallel
	// against a single sequence.
	align.pack_query_slots( SO::SLOT_0 | SO::SLOT_4, oligo(FORWARD) );
	align.pack_query_slots( SO::SLOT_1 | SO::SLOT_5, oligo(FORWARD).complement() );

	align.pack_query_slots( SO::SLOT_2 | SO::SLOT_6, oligo(REVERSE) );
	align.pack_query_slots( SO::SLOT_3 | SO::SLOT_7, oligo(REVERSE).complement() );

	const size_t num_background_amplicon = background_amplicons.size();
	
	for(unsigned int i = 0;i < num_background_amplicon;i += 2){
		
		align.pack_target_slots(SO::SLOT_0 | SO::SLOT_1, m_keys[background_amplicons[i].f]);
		align.pack_target_slots(SO::SLOT_2 | SO::SLOT_3, m_keys[background_amplicons[i].r]);
			
		if( (i + 1) < num_background_amplicon){
			
			align.pack_target_slots(SO::SLOT_4 | SO::SLOT_5, m_keys[background_amplicons[i + 1].f]);
			align.pack_target_slots(SO::SLOT_6 | SO::SLOT_7, m_keys[background_amplicons[i + 1].r]);
		}		

		align.align();

		/////////////////////////////////////////////////////////////////////////////////////////
		// Target i
		/////////////////////////////////////////////////////////////////////////////////////////
		float FpRm = align.score(0)*align.score(3)*f_norm*r_norm;
		float RpFm = align.score(1)*align.score(2)*f_norm*r_norm;

		if(m_opt.use_taq_mama){

			FpRm *= taq_mama_correction( Fp_taq_mama, align.target_last_two_aligned(0) )*
				taq_mama_correction( Rm_taq_mama, align.target_last_two_aligned(3) );


			RpFm *= taq_mama_correction( Rp_taq_mama, align.target_last_two_aligned(2) )*
				taq_mama_correction( Fm_taq_mama, align.target_last_two_aligned(1) );
		}

		float score = 0.0f;

		pair<int, int> plus_range;
		pair<int, int> minus_range;

		if(FpRm > RpFm){ // F(+) & R(-)

			score = sqrt(FpRm);

			plus_range = align.alignment_range_target(0);
			minus_range = align.alignment_range_target(3);

		}
		else{ // R(+) & F(-)

			score = sqrt(RpFm);

			plus_range = align.alignment_range_target(2);
			minus_range = align.alignment_range_target(1);
		}

		if(score >= m_opt.background_threshold){

			// This is a valid match
			m_match[ background_amplicons[i].index ] = true;
		}

		if( (i + 1) >= num_seq){
			continue;
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// Target i + 1
		/////////////////////////////////////////////////////////////////////////////////////////
		FpRm = align.score(4)*align.score(7)*f_norm*r_norm;
		RpFm = align.score(5)*align.score(6)*f_norm*r_norm;

		if(m_opt.use_taq_mama){

			FpRm *= taq_mama_correction( Fp_taq_mama, align.target_last_two_aligned(4) )*
				taq_mama_correction( Rm_taq_mama, align.target_last_two_aligned(7) );


			RpFm *= taq_mama_correction( Rp_taq_mama, align.target_last_two_aligned(6) )*
				taq_mama_correction( Fm_taq_mama, align.target_last_two_aligned(5) );
		}

		score = 0.0f;

		if(FpRm > RpFm){ // F(+) & R(-)

			score = sqrt(FpRm);

			plus_range = align.alignment_range_target(4);
			minus_range = align.alignment_range_target(7);

		}
		else{ // R(+) & F(-)

			score = sqrt(RpFm);

			plus_range = align.alignment_range_target(6);
			minus_range = align.alignment_range_target(5);
		}

		if(score >= m_opt.background_threshold){

			// This is a valid match
			m_match[ background_amplicons[i + 1].index ] = true;
		}
	}
}

void PCR::find_multiplex_background_match(BitSet &m_match, 
	const deque<Sequence> &m_seq, const Options &m_opt, std::ostream &m_out) const
{
	const unsigned int num_seq = m_seq.size();
	const float perfect_match_score = PERFECT_MATCH_SCORE;
	
	if(m_match.size() != num_seq){
		throw __FILE__ ":PCR::find_multiplex_background_match: |match| != num_seq";
	}

	float f_norm = perfect_match_score*oligo(FORWARD).size();
	float r_norm = perfect_match_score*oligo(REVERSE).size();
	
	if(f_norm > 0.0f){
		f_norm = 1.0f/f_norm;
	}

	if(r_norm > 0.0f){
		r_norm = 1.0f/r_norm;
	}
	
	pair<unsigned char, unsigned char> Fp_taq_mama;
	pair<unsigned char, unsigned char> Fm_taq_mama;
	pair<unsigned char, unsigned char> Rp_taq_mama;
	pair<unsigned char, unsigned char> Rm_taq_mama;
	
	if(m_opt.use_taq_mama){
		
		Fp_taq_mama = oligo(FORWARD).get_last_two();
		Fm_taq_mama = oligo(FORWARD).complement().get_last_two();
		Rp_taq_mama = oligo(REVERSE).get_last_two();
		Rm_taq_mama = oligo(REVERSE).complement().get_last_two();
	}
	
	SO::SeqOverlap align(SO::SeqOverlap::SmithWaterman, true /* is nucleic acid */);

	// Packing scheme:
	// slot 0: F   + seq i
	// slot 1: (F) + seq i
	// slot 2: R   + seq i
	// slot 3: (R) + seq i
	// slot 4: F   + seq i + 1
	// slot 5: (F) + seq i + 1
	// slot 6: R   + seq i + 1
	// slot 7: (R) + seq i + 1

	// Simple method #1. Search forward and reverse sequences in parallel
	// against a single sequence.
	align.pack_query_slots( SO::SLOT_0 | SO::SLOT_4, oligo(FORWARD) );
	align.pack_query_slots( SO::SLOT_1 | SO::SLOT_5, oligo(FORWARD).complement() );

	align.pack_query_slots( SO::SLOT_2 | SO::SLOT_6, oligo(REVERSE) );
	align.pack_query_slots( SO::SLOT_3 | SO::SLOT_7, oligo(REVERSE).complement() );

	for(unsigned int i = 0;i < num_seq;i += 2){

		if( (i + 1) < num_seq){

			align.pack_target_slots(SO::SLOT_0 | SO::SLOT_1 | SO::SLOT_2 | SO::SLOT_3, m_seq[i]);
			align.pack_target_slots(SO::SLOT_4 | SO::SLOT_5 | SO::SLOT_6 | SO::SLOT_7, m_seq[i + 1]);
		}
		else{
			align.pack_target(m_seq[i]);
		}

		align.align();

		/////////////////////////////////////////////////////////////////////////////////////////
		// Target i
		/////////////////////////////////////////////////////////////////////////////////////////
		float Fp_score = align.score(0)*f_norm;
		float Fm_score = align.score(1)*f_norm;
		
		float Rp_score = align.score(2)*r_norm;
		float Rm_score = align.score(3)*r_norm;
		
		if(m_opt.use_taq_mama){
			
			Fp_score *= taq_mama_correction( Fp_taq_mama, align.target_last_two_aligned(0) );
			Fm_score *= taq_mama_correction( Fm_taq_mama, align.target_last_two_aligned(1) );
			
			Rp_score *= taq_mama_correction( Rp_taq_mama, align.target_last_two_aligned(2) );
			Rm_score *= taq_mama_correction( Rm_taq_mama, align.target_last_two_aligned(3) );
		}

		// We have a match if *any* primer binds to a multiplex background amplicon
		if( (Fp_score >= m_opt.background_threshold) ||
		    (Fm_score >= m_opt.background_threshold) ||
		    (Rp_score >= m_opt.background_threshold) ||
		    (Rm_score >= m_opt.background_threshold) ){
		
			// This is a valid match
			m_match[i] = true;
		}

		if( (i + 1) >= num_seq){
			continue;
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// Target i + 1
		/////////////////////////////////////////////////////////////////////////////////////////
		Fp_score = align.score(4)*f_norm;
		Fm_score = align.score(5)*f_norm;
		
		Rp_score = align.score(6)*r_norm;
		Rm_score = align.score(7)*r_norm;

		if(m_opt.use_taq_mama){

			Fp_score *= taq_mama_correction( Fp_taq_mama, align.target_last_two_aligned(4) );
			Fm_score *= taq_mama_correction( Fm_taq_mama, align.target_last_two_aligned(5) );
			
			Rp_score *= taq_mama_correction( Rp_taq_mama, align.target_last_two_aligned(6) );
			Rm_score *= taq_mama_correction( Rm_taq_mama, align.target_last_two_aligned(7) );
		}

		// We have a match if *any* primer binds to a multiplex background amplicon
		if( (Fp_score >= m_opt.background_threshold) ||
		    (Fm_score >= m_opt.background_threshold) ||
		    (Rp_score >= m_opt.background_threshold) ||
		    (Rm_score >= m_opt.background_threshold) ){
		
			// This is a valid match
			m_match[i + 1] = true;
		}
	}
}
