#ifndef __PCRAMP
#define __PCRAMP

#include <string>
#include <deque>
#include <set>
#include "sequence.h"
#include "mpi_util.h"
#include "sort.h"

#define	PCRAMP_MAJOR_VERSION	"0"
#define	PCRAMP_MINOR_VERSION	"3"

#define	DEFAULT_MIN_TARGET_AMPLICON	80
#define	DEFAULT_MAX_TARGET_AMPLICON	200

#define	DEFAULT_MIN_BACKGROUND_AMPLICON	0
#define	DEFAULT_MAX_BACKGROUND_AMPLICON	2000

#define	DEFAULT_MIN_PRIMER	18
#define	DEFAULT_MAX_PRIMER	25

#define DEFAULT_MIN_PRIMER_TM	50.0f
#define DEFAULT_MAX_PRIMER_TM	75.0f
#define	DEFAULT_PRIMER_STRAND	900.0e-9f

#define	DEFAULT_SALT		0.05f
#define	DEFAULT_MAX_HAIRPIN	40.0f
#define	DEFAULT_MAX_DIMER	40.0f

#define	DEFAULT_DEGEN		1
#define	DEFAULT_NUM_TRIAL	1000
#define	DEFAULT_NUM_ASSAY	100

#define	DEFAULT_TARGET_WEIGHT		1.0f

#define DEFAULT_BACKGROUND_WEIGHT	1.0f

#define	DEFAULT_TARGET_THRESHOLD	1.0f
#define	DEFAULT_BACKGROUND_THRESHOLD	0.8f

#define	DEFAULT_PACK_MAX_DEGEN	256
#define	DEFAULT_PACK_MAX_GC	1.0f // Disabled
#define	DEFAULT_PACK_MIN_GC	0.0f // Disabled

// These values are now floating point values (since per-sequence
// weights are continuous values > 0.0
#define	DEFAULT_MIN_TARGET_COVER	0.0f
#define	DEFAULT_MAX_BACKGROUND_COVER	0.0f

#define	DEFAULT_SEARCH_THRESHOLD_MULTIPLIER	0.9

// When extracting the non-primer amplicon sequences for multiplex compatibility testing,
// include AMPLICON_PADDING extra 5' and 3' bases in the amplicon sequence (i.e. part of the
// primer binding site) to avoid unintended, weak primer binding that could produce competing
// parital amplicons
#define	MULTIPLEX_AMPLICON_PADDING	4

// For unix-like file systems
#define	PATH_SEPARATOR	'/'

// We need to protect any commas that appear in template variables
// if we are going to combine them with X Macros
#define SINGLE_ARG(...) __VA_ARGS__

struct Options
{	
	typedef enum {
		SILENT,
		VERBOSE,
		EVERYTHING,
		UNKNOWN_VERBOSITY
	} Verbosity;
	
	typedef enum {
		TEXT_OUTPUT,
		JSON_OUTPUT,
		UNKNOWN_OUTPUT
	} OutputFormat;
	
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized.
	#define OPTIONS_MEMBERS \
		VARIABLE(Verbosity, output_filter) \
		VARIABLE(OutputFormat, output_format) \
		VARIABLE(std::deque<std::string>, target_filename) \
		VARIABLE(std::deque<std::string>, background_filename) \
		VARIABLE(std::string, target_dir_prefix) \
		VARIABLE(std::string, background_dir_prefix) \
		VARIABLE(SINGLE_ARG(std::unordered_multimap<std::string, std::string>), target_groups) \
		VARIABLE(SINGLE_ARG(std::unordered_multimap<std::string, std::string>), background_groups) \
		VARIABLE(std::string, output_filename) \
		VARIABLE(unsigned int, degen) \
		VARIABLE(unsigned int, num_trial) \
		VARIABLE(unsigned int, num_assay) \
		VARIABLE(SINGLE_ARG(std::pair<int, int>), target_amplicon_range) \
		VARIABLE(SINGLE_ARG(std::pair<int, int>), background_amplicon_range) \
		VARIABLE(SINGLE_ARG(std::pair<int, int>), target_length_range) \
		VARIABLE(SINGLE_ARG(std::pair<int, int>), background_length_range) \
		VARIABLE(SINGLE_ARG(std::pair<int, int>), primer_range) \
		VARIABLE(SINGLE_ARG(std::pair<float, float>), primer_tm_range) \
		VARIABLE(float, max_hairpin) \
		VARIABLE(float, max_dimer) \
		VARIABLE(float, primer_strand) \
		VARIABLE(float, salt) \
		VARIABLE(float, target_weight) \
		VARIABLE(float, background_weight) \
		VARIABLE(float, target_search_multiplier) \
		VARIABLE(float, background_search_multiplier) \
		VARIABLE(float, target_threshold) \
		VARIABLE(float, background_threshold) \
		VARIABLE(float, min_target_cover) \
		VARIABLE(float, max_background_cover) \
		VARIABLE(unsigned int, seed) \
		VARIABLE(unsigned int, max_thread) \
		VARIABLE(unsigned int, pack_max_degen) \
		VARIABLE(float, pack_max_gc) \
		VARIABLE(float, pack_min_gc) \
		VARIABLE(std::deque<std::string>, target_ignore) \
		VARIABLE(std::deque<std::string>, background_ignore) \
		VARIABLE(bool, quit) \
		VARIABLE(bool, use_taq_mama) \
		VARIABLE(bool, top_down_search) \
		VARIABLE(bool, normalize_target_weight_per_file) \
		VARIABLE(bool, normalize_background_weight_per_file) \
		VARIABLE(bool, use_multiplex) \
		VARIABLE(bool, optimize_5) \
		VARIABLE(bool, optimize_3)
	
	#define VARIABLE(A, B) A B;
		OPTIONS_MEMBERS
	#undef VARIABLE
	
	Options();
	
	void load(int argc, char *argv[]);
	
	inline int min_oligo_length() const
	{
		return std::max(0, primer_range.first);
	};
	
	bool parse_json(const std::string &m_filename, const std::string &m_root_key,
		std::deque<std::string> &m_target_dir, std::string &m_target_dir_prefix,
		std::deque<std::string> &m_background_dir, std::string &m_background_dir_prefix);
	
	template<class T> friend size_t mpi_size(const T &m_obj);
	template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
		const T &m_obj);
	template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
		const T &m_obj);
};

template<> size_t mpi_size(const Options &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const Options &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, Options &m_obj);

struct Score
{
	// For multiplex assays, we will provide an additional reward for assays 
	// whose oligos are the *same* as the previously designed assays. The oligo_overlap
	// is the sum of the per-oligo fractional sequence similarities to the 
	// existing assay oligos.
	
	// Use X Macros (https://en.wikipedia.org/wiki/X_Macro) to 
	// ensure that structure variable are correctly serialized.
	#define SCORE_MEMBERS \
		VARIABLE(float, target_coverage) \
		VARIABLE(float, background_coverage) \
		VARIABLE(float, oligo_overlap)
	
	#define VARIABLE(A, B) A B;
		SCORE_MEMBERS
	#undef VARIABLE
	
	Score() :
		target_coverage(-1.0e6f), background_coverage(1.0e6f), oligo_overlap(0.0f)
	{
	};
	
	inline bool operator<(const Score &m_rhs) const
	{
		if( accuracy() == m_rhs.accuracy() ){
			return oligo_overlap < m_rhs.oligo_overlap;
		}
		
		return ( accuracy() < m_rhs.accuracy() );
	};
	
	inline bool operator>(const Score &m_rhs) const
	{
		if( accuracy() == m_rhs.accuracy() ){
			return oligo_overlap > m_rhs.oligo_overlap;
		}
		
		return ( accuracy() > m_rhs.accuracy() );
	};
	
	inline bool operator==(const Score &m_rhs) const
	{
		return ( ( accuracy() == m_rhs.accuracy() ) && 
		         (oligo_overlap == m_rhs.oligo_overlap) );
	};
	
	inline float accuracy() const
	{
		return target_coverage - background_coverage;
	};
	
	template<class T> friend size_t mpi_size(const T &m_obj);
	template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
		const T &m_obj);
	template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
		const T &m_obj);
};

template<> size_t mpi_size(const Score &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const Score &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, Score &m_obj);

template<class T>
inline void make_set(T &m_elements)
{
	// Use the gnu parallel sort 
	SORT( m_elements.begin(), m_elements.end() );

	// Return the unique elements
	m_elements.resize( std::unique( m_elements.begin(), m_elements.end() ) - m_elements.begin() );
};

template<class A, class B>
inline std::vector<A> keys(const MULTIMAP<A, B> &m_db)
{
	// Accumulate with a deque
	std::deque<A> ret;

	if( m_db.empty() ){
		return std::vector<A>();
	}
	
	ret.push_back(m_db.begin()->first);
	
	for(typename MULTIMAP<A, B>::const_iterator i = m_db.begin();i != m_db.end();++i){
	
		// Can't use x != y, so use !(x==y)
		if( !(ret.back() == i->first) ){
			ret.push_back(i->first);
		}
	}

	// Use the gnu parallel sort 
	SORT( ret.begin(), ret.end() );

	// Return the unique elements as a vector
	return std::vector<A>( ret.begin(),  std::unique( ret.begin(), ret.end() ) );
};

template<class A, class B>
inline std::vector<A> keys(const std::unordered_multimap<A, B> &m_db)
{
	// Accumulate with a deque
	std::deque<A> ret;

	if( m_db.empty() ){
		return std::vector<A>();
	}
	
	ret.push_back(m_db.begin()->first);
	
	for(typename std::unordered_multimap<A, B>::const_iterator i = m_db.begin();i != m_db.end();++i){
	
		// Can't use x != y, so use !(x==y)
		if( !(ret.back() == i->first) ){
			ret.push_back(i->first);
		}
	}

	// Use the gnu parallel sort 
	SORT( ret.begin(), ret.end() );

	// Return the unique elements as a vector
	return std::vector<A>( ret.begin(),  std::unique( ret.begin(), ret.end() ) );
};

// In filter.cpp
void filter_sequence(std::vector<Sequence> &m_seq, const Options &m_opt);

// In parse_fasta.cpp
void parse_fasta(const std::string &m_filename, std::deque<Sequence> &m_dbase, 
	const size_t &m_min_length, const size_t &m_max_length,
	const std::deque<std::string> &m_ignore);

void append_fasta_group(const std::string &m_filename, Sequence &m_seq,
	const size_t &m_min_length, const size_t &m_max_length,
	const size_t &m_num_pad, const std::deque<std::string> &m_ignore);

bool ignore_record(const std::string &m_defline, const std::deque<std::string> &m_ignore);

// In sample.cpp
int random_location(const int &m_start, const int &m_stop, unsigned int *m_seed_ptr);

// In options.cpp
std::string tolower(const std::string &m_str);

// background_support.cpp
std::vector<size_t> background_support(const std::deque<Sequence> &m_target_seq, 
	const std::deque<Sequence> &m_background_seq, const size_t &m_k);
	
#endif // __PCRAMP
