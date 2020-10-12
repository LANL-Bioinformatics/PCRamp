#include "assay.h"

using namespace std;

bool PCR::is_valid(const Oligo &m_oligo,const Word &m_trial_oligo, 
			NucCruc &m_melt, const Options &m_opt, 
			const bool &m_check_homo_dimer) const
{
	const double degen = m_trial_oligo.degeneracy();
	
	Word iter = m_trial_oligo.begin();
	
	m_melt.strand(m_opt.primer_strand/degen);
	
	do {
		const string seq = iter.str();
		
		float tm = m_melt.tm_pm_duplex(seq);

		if( (tm < m_opt.primer_tm_range.first) || (tm > m_opt.primer_tm_range.second) ){
			return false;
		}

		m_melt.set_query(seq);

		tm = m_melt.approximate_tm_hairpin();

		if(tm > m_opt.max_hairpin){
			return false;
		}


		if(m_check_homo_dimer){
		
			tm = m_melt.approximate_tm_homodimer();

			if(tm > m_opt.max_dimer){
				return false;
			}
		}
		
	} while( m_trial_oligo.next(iter) );
	
	return true;
}
