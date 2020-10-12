#ifndef __ASSAY
#define __ASSAY

#include "pcramp.h"
#include "bitset.h"
#include "nuc_cruc.h"
#include <ostream>
#include <deque>
#include <math.h>

// Enumerate the allowed oligos for PCR
enum {FORWARD, REVERSE, LAST_OLIGO};

// Sequence alignment scoring
#define	PERFECT_MATCH_SCORE	2.0f
#define MISMATCH_SCORE		-3.0f

// A score bonus to encourage oligo resuse when designing multiplex assays
#define	MULTIPLEX_OLIGO_REUSE_BONUS	10.0f

// The allowed optimization/relaxation moves
typedef enum {
	IncreaseDegeneracy, 
	DecreaseDegeneracy,
	Trim5, 
	Trim3, 
	Grow5, 
	Grow3,
	LastMove
} Move;
	
typedef unsigned char Oligo;

struct OligoMatch : public WordMatch
{
	unsigned int key_index;
	Oligo o;
	
	OligoMatch()
	{
	};
	
	OligoMatch(const Oligo &m_o, const unsigned int &m_w, const WordMatch &m_wm) :
		WordMatch(m_wm), key_index(m_w), o(m_o)
	{
	};
	
	inline bool operator<(const OligoMatch &m_rhs) const
	{
		// First sort by sequence id
		if(index < m_rhs.index){
			return true;
		}
		
		if(index > m_rhs.index){
			return false;
		}
		
		// Then sort by oligo location
		return (loc < m_rhs.loc);
	};
};

struct AmpliconBounds
{
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized.
	#define AMPLICON_BOUNDS_MEMBERS \
		VARIABLE(unsigned int, index) \
		VARIABLE(unsigned int, begin) \
		VARIABLE(unsigned int, end)
	
	#define VARIABLE(A, B) A B;
		AMPLICON_BOUNDS_MEMBERS
	#undef VARIABLE

	AmpliconBounds()
	{
		// Do nothing
	};

	AmpliconBounds(const unsigned int &m_index, const unsigned int &m_begin, const unsigned int &m_end):
		index(m_index), begin(m_begin), end(m_end)
	{
		if(begin > end){
			throw __FILE__ ":AmpliconBounds(): Amplicon begin > amplicon end";
		}
	};
};

template<> size_t mpi_size<AmpliconBounds>(const AmpliconBounds &m_obj);
template<> unsigned char* mpi_pack<AmpliconBounds>(unsigned char* m_ptr, const AmpliconBounds &m_obj);
template<> unsigned char* mpi_unpack<AmpliconBounds>(unsigned char* m_ptr, AmpliconBounds &m_obj);

// The following functions are in optimize.cpp
void match_words(std::deque<unsigned int> &m_match, const Word &m_seq, const std::vector<Word> &m_keys, const float &m_threshold);

void find_oligo_match(std::deque<OligoMatch> &m_match, const std::deque<unsigned int> &m_word_matches, const Oligo &m_oligo, 
	const Strand &m_strand, const std::vector<Word> &m_keys,
	const MULTIMAP<Word, WordMatch> &m_db,
	const std::deque<Sequence> &m_seq);

void update_identity(MAP<unsigned int, float> &m_ident, const Word &m_w, const std::vector<Word> &m_keys, bool m_use_taq_mama);

//////////////////////////////////////////////////////////////////////////////////////////////////////
// PCR assay
// |--F-->  <--R--|
//
// 1) Distance between 5' end of forward (F) and reverse (R) primers is the amplicon length, which 
// 	must be within allowed size range.
// 2) Forward and reverse primer oligos must not overlap.
// 3) Each oligo must be within the allow size range.
//////////////////////////////////////////////////////////////////////////////////////////////////////
class PCR
{
	private:
		Word f;
		Word r;
		
		struct PCROligos
		{
			// The sequence index
			unsigned int index;
			
			// The sequence weight
			float weight;
			
			// Indicies into the key array
			unsigned int f;
			unsigned int r;

			PCROligos()
			{
			};

			PCROligos(const unsigned int &m_index, const float &m_weight, 
				const unsigned int &m_f, const unsigned int &m_r) :
				index(m_index), weight(m_weight), f(m_f), r(m_r)
			{
			};

			inline bool operator<(const PCROligos &m_rhs) const
			{
				// First sort by sequence id
				if(index < m_rhs.index){
					return true;
				}

				if(index > m_rhs.index){
					return false;
				}

				// Then sort by oligo location
				if(f < m_rhs.f){
					return true;
				}

				if(f > m_rhs.f){
					return false;
				}

				return (r < m_rhs.r);
			};
		};
		
		struct sort_by_PCR_score
		{
			const MAP<unsigned int, float> &ref_f;
			const MAP<unsigned int, float> &ref_r;

			sort_by_PCR_score(const MAP<unsigned int, float> &m_f_identity, const MAP<unsigned int, float> &m_r_identity) :
				ref_f(m_f_identity), ref_r(m_r_identity)
			{
			};

			inline bool operator()(const PCROligos &m_lhs, const PCROligos &m_rhs) const
			{
				MAP<unsigned int, float>::const_iterator f_iter = ref_f.find(m_lhs.f);
				MAP<unsigned int, float>::const_iterator r_iter = ref_r.find(m_lhs.r);

				assert( ( f_iter != ref_f.end() ) && ( r_iter != ref_r.end() ) );
				
				const float score_lhs = sqrtf(f_iter->second * r_iter->second);
				
				f_iter = ref_f.find(m_rhs.f);
				r_iter = ref_r.find(m_rhs.r);

				assert( ( f_iter != ref_f.end() ) && ( r_iter != ref_r.end() ) );
				
				const float score_rhs = sqrtf(f_iter->second * r_iter->second);
				
				// Sort in descening order
				return (score_lhs > score_rhs);
			};
		};
		
		std::deque<PCROligos> target_amplicons;
		std::deque<PCROligos> background_amplicons;
		
		MAP<unsigned int, float> target_f_identity;
		MAP<unsigned int, float> target_r_identity;
		
		MAP<unsigned int, float> background_f_identity;
		MAP<unsigned int, float> background_r_identity;
		
		MAP<unsigned int, float> multiplex_background_f_identity;
		MAP<unsigned int, float> multiplex_background_r_identity;
		
		void collect_candidates(std::deque<PCROligos> &m_amplicons, 
				MAP<unsigned int, float> &m_f_identity, 
				MAP<unsigned int, float> &m_r_identity,
				const std::vector<Word> &m_keys, 
				const MULTIMAP<Word, WordMatch> &m_db, const std::deque<Sequence> &m_seq,
				const float &m_threshold,
				const std::pair<int, int> &m_amplicon_range) const;
		
		void sort_candidates(std::deque<PCROligos> &m_amplicons, 
				const MAP<unsigned int, float> &m_f_identity, 
				const MAP<unsigned int, float> &m_r_identity);
		
		void find_amplicon_match(std::deque<PCROligos> &m_amplicons, 
			std::deque<OligoMatch> &m_match, const Oligo &m_plus_oligo, const Oligo &m_minus_oligo, 
			const std::deque<Sequence> &m_seq, 
			const std::pair<int, int> &m_amplicon_range) const;
		
		void extract_amplicon_seq(std::deque<std::string> &m_amplicons, 
			std::deque<OligoMatch> &m_match, const Oligo &m_plus_oligo, const Oligo &m_minus_oligo, 
			const std::deque<Sequence> &m_seq, 
			const std::pair<int, int> &m_amplicon_range,
			std::deque<AmpliconBounds> *m_bounds_ptr = NULL);
			
		float compute_coverage(const std::deque<PCROligos> &m_amplicons, 
			const MAP<unsigned int, float> &m_f_identity, 
			const MAP<unsigned int, float> &m_r_identity, 
			const float &m_threshold) const;
		
		bool is_valid(const Oligo &m_oligo,const Word &m_trial_oligo, 
			NucCruc &m_melt, const Options &m_opt, const bool &m_check_homo_dimer) const;
	public:
		
		~PCR()
		{
			// Nothing to do!
		};
		
		inline void oligo(const Oligo &m_oligo, const Word &m_seq)
		{
			switch(m_oligo){
				case FORWARD:
					f = m_seq;
					break;
				case REVERSE:
					r = m_seq;
					break;
				default:
					throw __FILE__ ":PCR::oligo: Unknown oligo (1)";
			}
		};
		
		inline Word oligo(const Oligo &m_oligo) const
		{
			switch(m_oligo){
				case FORWARD:
					return f;
				case REVERSE:
					return r;
				default:
					throw __FILE__ ":PCR::oligo: Unknown oligo (2)";
			}
			
			// We should never get here
			return Word();
		};
		
		inline Oligo* assay_oligos() const
		{
			static Oligo ret[] = {FORWARD, REVERSE, LAST_OLIGO};
			
			return ret;
		};
		
		void copy_oligos(const PCR &m_rhs)
		{
			f = m_rhs.f;
			r = m_rhs.r;
		};
		
		inline void write(std::ostream &m_out) const
		{
			m_out << oligo(FORWARD) << '\t' <<  oligo(REVERSE)
				<< "\tD(F)=" << oligo(FORWARD).degeneracy() << ";D(R)=" 
				<< oligo(REVERSE).degeneracy();
		};

		inline void write_json(std::ostream &m_out) const
		{
			m_out << "\t\t\t\"forward primer\":{\n"
				<< "\t\t\t\t\"sequence\":\"" << oligo(FORWARD) << "\",\n"
				<< "\t\t\t\t\"degeneracy\":" << oligo(FORWARD).degeneracy() << "\n\t\t\t},\n"
				<< "\t\t\t\"reverse primer\":{\n"
				<< "\t\t\t\t\"sequence\":\"" << oligo(REVERSE) << "\",\n"
				<< "\t\t\t\t\"degeneracy\":" << oligo(REVERSE).degeneracy() << "\n\t\t\t},\n";
		};
		
		void write(std::ostream &m_out, const std::deque<PCR> &m_pool) const
		{
			// Since multiplexed assays are encouraged to reuse oligos, indicate
			// oligo reuse by writing the oligo in lower-case
			bool lower_f = false;
			bool lower_r = false;
			
			for(std::deque<PCR>::const_iterator i = m_pool.begin();i != m_pool.end();++i){
				
				if( std::max( oligo(FORWARD).max_overlap( i->oligo(FORWARD) ),
					      oligo(FORWARD).max_overlap( i->oligo(REVERSE) ) ) == 1.0f){
					lower_f = true;
				}
				
				if( std::max( oligo(REVERSE).max_overlap( i->oligo(FORWARD) ),
					      oligo(REVERSE).max_overlap( i->oligo(REVERSE) ) ) == 1.0f){
					lower_r = true;
				}
			}
			
			if(lower_f){			
				m_out << tolower( oligo(FORWARD).str() );
			}
			else{
				m_out << oligo(FORWARD);
			}
			
			m_out << '\t';
			
			if(lower_r){			
				m_out << tolower( oligo(REVERSE).str() );
			}
			else{
				m_out << oligo(REVERSE);
			}
			
			m_out << "\tD(F)=" << oligo(FORWARD).degeneracy() << ";D(R)=" 
				<< oligo(REVERSE).degeneracy();
		};
		
		void write_json(std::ostream &m_out, const std::deque<PCR> &m_pool) const
		{
			// Since multiplexed assays are encouraged to reuse oligos, indicate
			// oligo reuse by writing the oligo in lower-case
			bool lower_f = false;
			bool lower_r = false;
			
			for(std::deque<PCR>::const_iterator i = m_pool.begin();i != m_pool.end();++i){
				
				if( std::max( oligo(FORWARD).max_overlap( i->oligo(FORWARD) ),
					      oligo(FORWARD).max_overlap( i->oligo(REVERSE) ) ) == 1.0f){
					lower_f = true;
				}
				
				if( std::max( oligo(REVERSE).max_overlap( i->oligo(FORWARD) ),
					      oligo(REVERSE).max_overlap( i->oligo(REVERSE) ) ) == 1.0f){
					lower_r = true;
				}
			}
			
			m_out << "\t\t\t\"forward primer\":{\n"
				<< "\t\t\t\t\"sequence\":\"" << oligo(FORWARD) << "\",\n"
				<< "\t\t\t\t\"degeneracy\":" << oligo(FORWARD).degeneracy() << ",\n"
				<< "\t\t\t\t\"recycled\":" << (lower_f ? "True" : "False") << "\n"
				<< "\t\t\t},\n"
				<< "\t\t\t\"reverse primer\":{\n"
				<< "\t\t\t\t\"sequence\":\"" << oligo(REVERSE) << "\",\n"
				<< "\t\t\t\t\"degeneracy\":" << oligo(REVERSE).degeneracy() << ",\n"
				<< "\t\t\t\t\"recycled\":" << (lower_r ? "True" : "False") << "\n"
				<< "\t\t\t},\n";
		};
		
		inline bool operator<(const PCR &m_rhs) const
		{
			if(f < m_rhs.f){
				return true;
			}
			
			if(f > m_rhs.f){
				return false;
			}
			
			return (r < m_rhs.r);
		};
		
		inline bool operator==(const PCR &m_rhs) const
		{
			return (f == m_rhs.f) && (r == m_rhs.r);
		};
		
		inline void center()
		{
			f.center();
			r.center();
		};
		
		void collect_target_candidates(
			const std::vector<Word> &m_keys, 
			const MULTIMAP<Word, WordMatch> &m_db, const std::deque<Sequence> &m_seq, const Options &m_opt)
		{
			collect_candidates(target_amplicons, target_f_identity, target_r_identity,
				m_keys, m_db, m_seq, 
				m_opt.target_threshold*m_opt.target_search_multiplier,
				m_opt.target_amplicon_range);
		};
		
		void collect_background_candidates(
			const std::vector<Word> &m_keys, 
			const MULTIMAP<Word, WordMatch> &m_db, const std::deque<Sequence> &m_seq, const Options &m_opt)
		{
			if( !m_keys.empty() ){
				collect_candidates(background_amplicons, background_f_identity, background_r_identity,
					m_keys, m_db, m_seq, 
					m_opt.background_threshold*m_opt.background_search_multiplier,
					m_opt.background_amplicon_range);
			}
		};
		
		void collect_multiplex_background_candidates(
			const std::vector<Word> &m_keys, 
			const MULTIMAP<Word, WordMatch> &m_db, const std::deque<Sequence> &m_seq, const Options &m_opt);
		
		void sort_target_candidates()
		{
			sort_candidates(target_amplicons, target_f_identity, target_r_identity);
		};
		
		bool maximize_degeneracy(const std::vector<Word> &m_keys, NucCruc &m_melt, const Options &m_opt);
		
		float max_dimer_tm(NucCruc &m_melt, const Options &m_opt) const;
		
		inline void update_target_candidates(const std::vector<Word> &m_keys, const bool &m_use_taq_mama)
		{
			update_identity(target_f_identity, f, m_keys, m_use_taq_mama);
			update_identity(target_r_identity, r, m_keys, m_use_taq_mama);
		};
		
		inline void update_background_candidates(const std::vector<Word> &m_keys, const bool &m_use_taq_mama)
		{
			update_identity(background_f_identity, f, m_keys, m_use_taq_mama);
			update_identity(background_r_identity, r, m_keys, m_use_taq_mama);
		};
		
		inline void update_multiplex_background_candidates(const std::vector<Word> &m_keys, const bool &m_use_taq_mama)
		{
			update_identity(multiplex_background_f_identity, f, m_keys, m_use_taq_mama);
			update_identity(multiplex_background_r_identity, r, m_keys, m_use_taq_mama);
		};
		
		inline float compute_target_coverage(const float &m_threshold) const
		{
			return compute_coverage(target_amplicons, target_f_identity, target_r_identity, m_threshold);
		};
		
		inline float compute_background_coverage(const float &m_threshold) const
		{
			return compute_coverage(background_amplicons, background_f_identity, background_r_identity, m_threshold);
		};
		
		float compute_multiplex_background_coverage(const float &m_threshold) const;
		
		void find_target_match(BitSet &m_match, const std::vector<Word> &m_keys, 
			const MULTIMAP<Word, WordMatch> &m_db, const std::deque<Sequence> &m_seq, const Options &m_opt);
		
		void find_background_match(BitSet &m_match, const std::vector<Word> &m_keys, 
			const MULTIMAP<Word, WordMatch> &m_db, const std::deque<Sequence> &m_seq, 
			const Options &m_opt, std::ostream &m_out);
		
		void find_multiplex_background_match(BitSet &m_match, 
			const std::deque<Sequence> &m_seq, const Options &m_opt, std::ostream &m_out) const;
		
		bool multiplex_compatible(NucCruc &m_melt, 
			const Options &m_opt, const PCR &m_assay) const;

		std::pair<Word, Score> increase_degeneracy(const Oligo &m_oligo, 
				const std::vector<Word> &m_target_keys, const float &m_target_threshold,
				const std::vector<Word> &m_background_keys, const float &m_background_threshold,
				const std::vector<Word> &m_multiplex_background_keys,
				const Score &m_score_threshold, NucCruc &m_melt, 
				const std::deque<PCR> &m_pool, const Options &m_opt);
		
		std::pair<Word, Score> decrease_degeneracy(const Oligo &m_oligo, 
				const std::vector<Word> &m_target_keys, const float &m_target_threshold,
				const std::vector<Word> &m_background_keys, const float &m_background_threshold,
				const std::vector<Word> &m_multiplex_background_keys,
				const Score &m_score_threshold, NucCruc &m_melt, 
				const std::deque<PCR> &m_pool, const Options &m_opt);
		
		std::pair<Word, Score> trim5(const Oligo &m_oligo, 
				const std::vector<Word> &m_target_keys, const float &m_target_threshold,
				const std::vector<Word> &m_background_keys, const float &m_background_threshold,
				const std::vector<Word> &m_multiplex_background_keys,
				const Score &m_score_threshold, NucCruc &m_melt, 
				const std::deque<PCR> &m_pool, const Options &m_opt);
				
		std::pair<Word, Score> trim3(const Oligo &m_oligo, 
				const std::vector<Word> &m_target_keys, const float &m_target_threshold,
				const std::vector<Word> &m_background_keys, const float &m_background_threshold,
				const std::vector<Word> &m_multiplex_background_keys,
				const Score &m_score_threshold, NucCruc &m_melt, 
				const std::deque<PCR> &m_pool, const Options &m_opt);
		
		std::pair<Word, Score> grow5(const Oligo &m_oligo, 
				const std::vector<Word> &m_target_keys, const float &m_target_threshold,
				const std::vector<Word> &m_background_keys, const float &m_background_threshold,
				const std::vector<Word> &m_multiplex_background_keys,
				const Score &m_score_threshold, NucCruc &m_melt, 
				const std::deque<PCR> &m_pool, const Options &m_opt);
		
		std::pair<Word, Score> grow3(const Oligo &m_oligo, 
				const std::vector<Word> &m_target_keys, const float &m_target_threshold,
				const std::vector<Word> &m_background_keys, const float &m_background_threshold,
				const std::vector<Word> &m_multiplex_background_keys,
				const Score &m_score_threshold, NucCruc &m_melt, 
				const std::deque<PCR> &m_pool, const Options &m_opt);
			
		void random_assay(const std::deque<Sequence> &m_seq, 
			NucCruc &m_melt, const Options &m_opt, unsigned int &m_seed,
			std::ostream &m_out);
		
		inline std::string packed_string() const
		{
			std::string ret;
			
			f.write(ret);
			ret.push_back('|'); // Delimiter to make packed string uniquely map to assays
			r.write(ret);
			
			return ret;
		};
		
		inline double total_degeneracy() const
		{
			return f.degeneracy() + r.degeneracy();
		};
		
		float compute_oligo_overlap(const std::deque<PCR> &m_pool) const;
				
		std::deque<Sequence> collect_unique_amplicons(
			const std::vector<Word> &m_keys, 
			const MULTIMAP<Word, WordMatch> &m_db, 
			const std::deque<Sequence> &m_seq, 
			const float &m_threshold,
			const std::pair<int, int> &m_amplicon_range,
			std::deque<AmpliconBounds> *m_bounds_ptr = NULL);
		
		inline size_t mpi_size() const
		{
			return f.mpi_size() + r.mpi_size();
		};
		
		inline unsigned char* mpi_pack(unsigned char* m_ptr) const
		{
			m_ptr = f.mpi_pack(m_ptr);
			m_ptr = r.mpi_pack(m_ptr);
			
			return m_ptr;
		
		};
		
		inline unsigned char* mpi_unpack(unsigned char* m_ptr)
		{
			m_ptr = f.mpi_unpack(m_ptr);
			m_ptr = r.mpi_unpack(m_ptr);
			
			return m_ptr;
		};
};

// In optimize.cpp
Score optimize(PCR &m_assay, const std::vector<Move> &m_moves,
		const std::vector<Word> &m_target_keys, 
		const MULTIMAP<Word, WordMatch> &m_target_db, 
		const std::deque<Sequence> &m_target_seq, 
		const std::vector<Word> &m_background_keys, 
		const MULTIMAP<Word, WordMatch> &m_background_db, 
		const std::deque<Sequence> &m_background_seq,
		const std::vector<Word> &m_multiplex_background_keys, 
		const MULTIMAP<Word, WordMatch> &m_multiplex_background_db, 
		const std::deque<Sequence> &m_multiplex_background_seq, 
		const std::deque<PCR> &m_pool,
		const Options &m_opt,
		std::ostream &m_out);

std::pair<Word, Score> optimization_move(const Move &m_move,
	const Oligo &m_oligo, 
	PCR &m_assay, 
	const std::vector<Word> &m_target_keys, const float &m_target_threshold,
	const std::vector<Word> &m_background_keys, const float &m_background_threshold,
	const std::vector<Word> &m_multiplex_background_keys,
	const Score &m_score_threshold, NucCruc &m_melt, 
	const std::deque<PCR> &m_pool, 
	const Options &m_opt);

bool make_degenerate(PCR &m_assay, 
		const std::vector<Word> &m_target_keys, 
		const MULTIMAP<Word, WordMatch> &m_target_db, 
		const std::deque<Sequence> &m_target_seq, 
		NucCruc &m_melt, const Options &m_opt,
		std::ostream &m_out);

// In select_words.cpp
void select_words(MULTIMAP<Word, WordMatch> &m_dst, 
	const MULTIMAP<Word, WordMatch> &m_src,
	const std::vector<PCR> &m_candidates, 
	const bool &m_optimize_5, const bool &m_optimize_3, 
	const float &m_threshold);
	
#endif // __ASSAY
