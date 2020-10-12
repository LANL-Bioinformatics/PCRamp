#include "pcramp.h"
#include <limits.h>
#include <stdlib.h>
#include <getopt.h>

#include <mpi.h>

#include <iostream>
#include <fstream>

// For testing and navigating directories
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>

#include "JSON.h"

using namespace std;

Options::Verbosity parse_verbosity(string m_buffer);

// Enumerate the allowed fasta file extensions
const char* allowed_fasta_extions [] = {
	".fna", ".fna.gz",
	".fasta", ".fasta.gz",
	".fa", ".fa.gz",
	NULL
};

// Global MPI variables
extern int numtasks;
extern int id;

bool find_groups(const string &m_path, unordered_multimap<string, string> &m_group);
bool find_file_extension(const string &m_path, const char** m_ext);
void strip_trailing_path_separtor(deque<string> &m_str);
deque<string> parse_keys(const string &m_key_str);
string replace_special_with(char m_replacement, const string &m_str);

Options::Options()
{
	// Set the default options
	degen = DEFAULT_DEGEN;
	num_trial = DEFAULT_NUM_TRIAL;
	num_assay = DEFAULT_NUM_ASSAY;
	output_filter = VERBOSE;
	
	output_format = TEXT_OUTPUT;
	
	target_amplicon_range = make_pair(DEFAULT_MIN_TARGET_AMPLICON, DEFAULT_MAX_TARGET_AMPLICON);
	background_amplicon_range = make_pair(DEFAULT_MIN_BACKGROUND_AMPLICON, DEFAULT_MAX_BACKGROUND_AMPLICON);
	
	target_length_range = make_pair(0, INT_MAX);
	background_length_range = make_pair(0, INT_MAX);
	
	primer_range = make_pair(DEFAULT_MIN_PRIMER, DEFAULT_MAX_PRIMER);
	
	primer_tm_range = make_pair(DEFAULT_MIN_PRIMER_TM, DEFAULT_MAX_PRIMER_TM);
	
	primer_strand = DEFAULT_PRIMER_STRAND;

	salt = DEFAULT_SALT;
	max_hairpin = DEFAULT_MAX_HAIRPIN;
	max_dimer = DEFAULT_MAX_DIMER;
	
	quit = false;
	use_taq_mama = false;
	top_down_search = false;
	normalize_target_weight_per_file = false;
	normalize_background_weight_per_file = false;
	use_multiplex = true;
	optimize_5 = false;
	optimize_3 = false;
	
	target_weight = DEFAULT_TARGET_WEIGHT;
	background_weight = DEFAULT_BACKGROUND_WEIGHT;
	
	target_search_multiplier = DEFAULT_SEARCH_THRESHOLD_MULTIPLIER;
	background_search_multiplier = DEFAULT_SEARCH_THRESHOLD_MULTIPLIER;
	target_threshold = DEFAULT_TARGET_THRESHOLD;
	background_threshold = DEFAULT_BACKGROUND_THRESHOLD;
	min_target_cover = DEFAULT_MIN_TARGET_COVER;
	max_background_cover = DEFAULT_MAX_BACKGROUND_COVER;
		
	seed = 0; // A seed of 0 triggers a time-based seed
	
	// A max_thread of 0 uses all available threads
	max_thread = 0;
	
	pack_max_degen = DEFAULT_PACK_MAX_DEGEN;
	pack_max_gc = DEFAULT_PACK_MAX_GC;
	pack_min_gc = DEFAULT_PACK_MIN_GC;
}

void Options::load(int argc, char *argv[])
{
	
	bool print_usage = (argc <= 1);
	
	// Command line options:
	// -t <fasta file of individual targets>
	// -T <root directory of target subdirectories>
	// [--T.prefix <directory prefix for target genomes>]
	// [-b <fasta file of individual backgrounds>]
	// [-B <root directory of background subdirectories>]
	// [--B.prefix <directory prefix for background genomes>]
	// -o <output file>
	// [-d <max degeneracy> (default is 1)]
	// [--count <total number of PCR assays to produce>]
	// [--optimize.top-down Search using maximally degenerate inital primers (default is bottom up)]
	// [--salt <salt concentration>]
	// [--primer.hairpin <max oligo hairpin  melting temperature>]
	// [--primer.dimer <max oligo hetero/homo-dimer  melting temperature>]
	// [--optimize.5 (enable 5' oligo search to optimize primer coverage)]
	// [--no-optimize.5 (disable 5' oligo search to optimize primer coverage)]
	// [--optimize.3 (enable 3' oligo search to optimize primer coverage)]
	// [--no-optimize.3 (disable 3' oligo search to optimize primer coverage)]
	// [--target.amplicon.min <minimum amplicon length> (default is 80)]
	// [--target.amplicon.max <maximum amplicon length> (default is 200)]
	// [--target.weight <target weight> (default is 1.0)]
	// [--target.threshold <target detection threshold> (default is 1.0)]
	// [--target.cover <minimum per primer pair coverage> (default is 1)]
	// [--target.ignore <defline key word to exclude a sequence>]
	// [--target.normalize (normalize weights per fasta file)]
	// [--target.size.min <length> (minimum input target length in bp)]
	// [--target.size.max <length> (maximum input target length in bp)]
	// [--target.prefix <directory prefix for target genomes>] (synonym for --T.prefix)
	// [--background.amplicon.min <minimum amplicon length> (default is 80)]
	// [--background.amplicon.max <maximum amplicon length> (default is 200)]
	// [--background.weight <background weight> (default is 1.0)]
	// [--background.threshold <background detection threshold> (default is 1.0)]
	// [--background.cover <maximum per primer pair coverage> (default is 0)]
	// [--background.ignore <defline key word to exclude a sequence>]
	// [--background.normalize (normalize weights per fasta file)]
	// [--background.size.min <length> (minimum input background length in bp)]
	// [--background.size.max <length> (maximum input background length in bp)]
	// [--background.fast (Only use a subset of background sequences)]
	// [--background.prefix <directory prefix for background genomes>] (synonym for --B.prefix)
	// [--primer.size.min <minimum primer length> (default is 18)]
	// [--primer.size.max <maximum primer length> (default is 25)]
	// [--primer.tm.min <minimum primer melting temperature>]
	// [--primer.tm.max <maximum primer melting temperature>]
	// [--primer.strand <primer concentration>]
	// [--primer.taq-mama (use Taq MAMA rules for terminal primer mismatches)]
	// [--trial <number of trials> (default is 1000)]
	// [--seed <random number seed> (default is time-based)]
	// [--thread <maximum number of OpenMP threads> (default is all)]
	// [-v <verbosity level: silent, verbose, everything>]
	// [--pack.degen.max <max degen when packing words> (default is 256)]
	// [--pack.gc.max <max fractional GC content when packing words> (default is 1.0, disabled)]
	// [--pack.gc.min <min fractional GC content when packing words> (default is 0.0, disabled)]
	// [--input.prefix <directory prefix for both target and background input genomes>]
	// [--json <JSON file of configuration values>]
	// [--json.root <root key for JSON configuration file> ("xxx|yyy|zzz|...")]
	
	const char* options = "t:T:b:B:o:d:v:?h";
	int config_opt = 0;
	int long_index = 0;

	struct option long_opts[] = {
		{"target.amplicon.min", true, &config_opt, 2},
		{"target.amplicon.max", true, &config_opt, 3},
		{"primer.size.min", true, &config_opt, 4},
		{"primer.size.max", true, &config_opt, 5},
		{"background.amplicon.min", true, &config_opt, 8},
		{"background.amplicon.max", true, &config_opt, 9},
		{"seed", true, &config_opt, 10},
		{"target.weight", true, &config_opt, 11},
		{"background.weight", true, &config_opt, 12},
		{"target.threshold", true, &config_opt, 13},
		{"background.threshold", true, &config_opt, 14},
		{"target.cover", true, &config_opt, 15},
		{"background.cover", true, &config_opt, 16},
		{"primer.taq-mama", false, &config_opt, 17},
		{"trial", true, &config_opt, 18},
		{"pack.degen.max", true, &config_opt, 19},
		{"pack.gc.min", true, &config_opt, 20},
		{"pack.gc.max", true, &config_opt, 21},
		{"count", true, &config_opt, 23},
		{"target.search", true, &config_opt, 24},
		{"optimize.top-down", false, &config_opt, 26},
		{"target.ignore", true, &config_opt, 27},
		{"background.ignore", true, &config_opt, 28},
		{"primer.tm.min", true, &config_opt, 29},
		{"primer.tm.max", true, &config_opt, 30},
		{"primer.strand", true, &config_opt, 33},
		{"salt", true, &config_opt, 35},
		{"primer.hairpin", true, &config_opt, 36},
		{"primer.dimer", true, &config_opt, 37},
		{"target.normalize", false, &config_opt, 39},
		{"background.normalize", false, &config_opt, 40},
		{"thread", true, &config_opt, 41},
		//{"multiplex", false, &config_opt, 42}, <-- default is true
		{"target.size.min", true, &config_opt, 43},
		{"target.size.max", true, &config_opt, 44},
		{"background.size.min", true, &config_opt, 45},
		{"background.size.max", true, &config_opt, 46},
		{"background.search", true, &config_opt, 47},
		{"optimize.5", false, &config_opt, 48},
		{"no-optimize.5", false, &config_opt, 49},
		{"optimize.3", false, &config_opt, 50},
		{"no-optimize.3", false, &config_opt, 51},
		{"json", true, &config_opt, 53},
		{"json.root", true, &config_opt, 54},
		{"target.prefix", true, &config_opt, 55},
		{"T.prefix", true, &config_opt, 55}, // Synonym for target.prefix
		{"background.prefix", true, &config_opt, 56},
		{"B.prefix", true, &config_opt, 56}, // Synonym for background.prefix
		{"input.prefix", true, &config_opt, 57},
		{"o.text", false, &config_opt, 58},
		{"o.json", false, &config_opt, 59},
		{0,0,0,0} // Terminate options list
	};

	int opt_code;
	opterr = 0;

	deque<string> target_dir;
	deque<string> background_dir;
	
	string json_file;
	string json_root;
	
	string dir_prefix;
	
	while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){

		switch( opt_code ){
			case 0:
				
				if(config_opt == 2){ // target.amplicon.min

					target_amplicon_range.first = atoi(optarg);
					
					if(target_amplicon_range.first < 0){
						
						cerr << "Please specify a target.amplicon.min >= 0" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}

				if(config_opt == 3){ // target.amplicon.max

					target_amplicon_range.second = atoi(optarg);
					
					if(target_amplicon_range.second < 0){
						
						cerr << "Please specify a target.amplicon.max >= 0" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}

				if(config_opt == 4){ // primer.size.min

					primer_range.first = atoi(optarg);
					
					if(primer_range.first < 0){
						
						cerr << "Please specify a primer.min >= 0" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}

				if(config_opt == 5){ // primer.size.max

					primer_range.second = atoi(optarg);
					
					if(primer_range.second < 0){
						
						cerr << "Please specify a primer.max >= 0" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}
								
				if(config_opt == 8){ // background.amplicon.min

					background_amplicon_range.first = atoi(optarg);
					
					if(background_amplicon_range.first < 0){
						
						cerr << "Please specify a background.amplicon.min >= 0" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}

				if(config_opt == 9){ // background.amplicon.max

					background_amplicon_range.second = atoi(optarg);
					
					if(background_amplicon_range.second < 0){
						
						cerr << "Please specify a background.amplicon.max >= 0" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}
				
				if(config_opt == 10){ // seed

					seed = abs( atoi(optarg) );
					break;
				}
				
				if(config_opt == 11){ // target.weight

					target_weight = atof(optarg);
					
					if(target_weight < 1.0){
					
						cerr << "Please specify a valid target.weight value (>= 1.0)" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}
				
				if(config_opt == 12){ // background.weight

					background_weight = atof(optarg);
					
					if(background_weight < 0.0f){
					
						cerr << "Please specify a valid background.weight value (>= 0.0)" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}
				
				if(config_opt == 13){ // target.threshold

					target_threshold = atof(optarg);
					
					if( (target_threshold < 0.0f) || (target_threshold > 1.0f) ){
					
						cerr << "Please specify a valid target.threshold value (0 <= PCR <= 1)" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}
				
				if(config_opt == 14){ // background.threshold

					background_threshold = atof(optarg);
					
					if( (background_threshold < 0.0f) || (background_threshold > 1.0f) ){
					
						cerr << "Please specify a valid background.threshold value (0 <= PCR <= 1)" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}
				
				if(config_opt == 15){ // target.cover

					min_target_cover = atof(optarg);
					
					if(min_target_cover < 0.0f){
					
						cerr << "Please specify a valid target.cover value (>= 0)" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}
				
				if(config_opt == 16){ // background.cover

					max_background_cover = atof(optarg);
					
					if(max_background_cover < 0.0f){
					
						cerr << "Please specify a valid background.cover value (>= 0)" << endl;
						
						quit = true;
						return;
					}
										
					break;
				}
				
				if(config_opt == 17){ // primer.taq-mama

					use_taq_mama = true;										
					break;
				}
				
				if(config_opt == 18){ // trial

					// Read the number of trials as a float to allow the use
					// of scientific notation when specifying large values for
					// the number of trials.
					num_trial = (unsigned int)( fabs( atof(optarg) ) );
					
					if(num_trial == 0){
					
						cerr << "Please enter --trial > 0" << endl;
						
						quit = true;
						return;
					}
														
					break;
				}
				
				if(config_opt == 19){ // pack.degen.max

					pack_max_degen = abs( atoi(optarg) );
					
					if(pack_max_degen == 0){
					
						cerr << "Please enter --pack.max_degen > 0" << endl;
						
						quit = true;
						return;
					}
														
					break;
				}
				
				if(config_opt == 20){ // pack.gc.min

					pack_min_gc = atof(optarg);
					
					if( (pack_min_gc < 0.0f) || (pack_min_gc > 1.0f) ){
						
						cerr << "Please enter --pack.min_gc >= 0 and <= 1.0" << endl;
						
						quit = true;
						return;
					}
														
					break;
				}
				
				if(config_opt == 21){ // pack.gc.max

					pack_max_gc = atof(optarg);
					
					if( (pack_max_gc < 0.0f) || (pack_max_gc > 1.0f) ){
						
						cerr << "Please enter --pack.max_gc >= 0 and <= 1.0" << endl;
						
						quit = true;
						return;
					}
														
					break;
				}
								
				if(config_opt == 23){ // count

					num_assay = abs( atoi(optarg) );
					
					if( num_assay == 0 ){
						
						cerr << "Please enter --count >= 1" << endl;
						
						quit = true;
						return;
					}
														
					break;
				}
				
				if(config_opt == 24){ // target.search

					target_search_multiplier = atof(optarg);
					
					if( (target_search_multiplier <= 0.0) || (target_search_multiplier > 1.0) ){
						
						cerr << "Please enter 0 < target.search <= 1" << endl;
						
						quit = true;
						return;
					}
														
					break;
				}
				
				if(config_opt == 26){ // optimize.top-down

					top_down_search = true;						
					break;
				}
				
				if(config_opt == 27){ // target.ignore

					target_ignore.push_back( tolower(optarg) );
					break;
				}
				
				if(config_opt == 28){ // background.ignore

					background_ignore.push_back( tolower(optarg) );
					break;
				}
				
				if(config_opt == 29){ // primer.tm.min
				
					primer_tm_range.first = atof(optarg);
					break;
				}
				
				if(config_opt == 30){ // primer.tm.max
				
					primer_tm_range.second = atof(optarg);
					break;
				}
								
				if(config_opt == 33){ // primer.strand
				
					primer_strand = atof(optarg);
					break;
				}
								
				if(config_opt == 35){ // salt
				
					salt = atof(optarg);
					break;
				}
				
				if(config_opt == 36){ // primer.hairpin
				
					max_hairpin = atof(optarg);
					break;
				}
				
				if(config_opt == 37){ // primer.dimer
				
					max_dimer = atof(optarg);
					break;
				}
								
				if(config_opt == 39){ // target.normalize
				
					normalize_target_weight_per_file = true;
					break;
				}
				
				if(config_opt == 40){ // background.normalize
				
					normalize_background_weight_per_file = true;
					break;
				}
				
				if(config_opt == 41){ // thread
				
					max_thread = abs( atoi(optarg) );
					break;
				}
				
				// Multiplex assay design is now the default
				//if(config_opt == 42){ // multiplex
				//
				//	use_multiplex = true;
				//	break;
				//}
				
				if(config_opt == 43){ // target.size.min
				
					target_length_range.first = atoi(optarg);
					break;
				}
				
				if(config_opt == 44){ // target.size.max
				
					target_length_range.second = atoi(optarg);
					break;
				}
				
				if(config_opt == 45){ // background.size.min
				
					background_length_range.first = atoi(optarg);
					break;
				}
				
				if(config_opt == 46){ // background.size.max
				
					background_length_range.second = atoi(optarg);
					break;
				}
				
				if(config_opt == 47){ // background.search
				
					background_search_multiplier = atof(optarg);
					
					if( (background_search_multiplier <= 0.0) || (background_search_multiplier > 1.0) ){
						
						cerr << "Please enter 0 < background.search <= 1" << endl;
						
						quit = true;
						return;
					}
					
					break;
				}
				
				if(config_opt == 48){ // optimize.5
					
					optimize_5 = true;
					break;
				}
				
				if(config_opt == 49){ // no-optimize.5
					
					optimize_5 = false;
					break;
				}
				
				if(config_opt == 50){ // optimize.3
					
					optimize_3 = true;
					break;
				}
				
				if(config_opt == 51){ // no-optimize.3
					
					optimize_3 = false;
					break;
				}
								
				if(config_opt == 53){ // json
					
					json_file = optarg;
					break;
				}
				
				if(config_opt == 54){ // json.root
					
					json_root = optarg;
					break;
				}
				
				if(config_opt == 55){ // target.prefix or T.prefix
					
					target_dir_prefix = optarg;
					break;
				}
				
				if(config_opt == 56){ // background.prefix or B.prefix
					
					background_dir_prefix = optarg;
					break;
				}
				
				if(config_opt == 57){ // input.prefix
					
					dir_prefix = optarg;
					break;
				}
				
				if(config_opt == 58){ // o.text
					
					output_format = TEXT_OUTPUT;
					break;
				}
				
				if(config_opt == 59){ // o.json
					
					output_format = JSON_OUTPUT;
					break;
				}
				
				cerr << "Unknown command line flag!" << endl;
				
				quit = true;
				return;
				
			case 't':
				target_filename.push_back(optarg);
				break;
			case 'T':
				target_dir.push_back(optarg);
				break;
			case 'b':
				background_filename.push_back(optarg);
				break;
			case 'B':
				background_dir.push_back(optarg);
				break;
			case 'o':
				output_filename = optarg;
				break;
			case 'd':
				degen = abs( atoi(optarg) );
				break;
			case 'v':
				output_filter = parse_verbosity(optarg);
				
				if(output_filter == UNKNOWN_VERBOSITY){
				
					cerr <<  "Please enter a valid verbosity flag: \"silent\", \"verbose\", \"everything\"" 
						<< endl;
						
					quit = true;
					return;
				}
				
				break;
			case 'h':
			case '?':
				print_usage = true;
				break;
			default:
				cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
				
				quit = true;
				return;
		};
	}
	
	if(print_usage){
		
		cerr << "PCRamp version " 
			<< PCRAMP_MAJOR_VERSION << "." 
			<< PCRAMP_MINOR_VERSION << endl;
		cerr << "Usage:" << endl;
		cerr << "\t-t <target fasta file>" << endl;
		cerr << "\t-T <root directory of target subdirectories>" << endl;
		cerr << "\t[--T.prefix <directory prefix for target genomes>]" << endl;
		cerr << "\t[-b <background fasta file>]" << endl;
		cerr << "\t[-B <root directory of background subdirectories>]" << endl;
		cerr << "\t[--B.prefix <directory prefix for background genomes>]" << endl;
		cerr << "\t-o <output file>" << endl;
		cerr << "\t[--o.text (output results as [poorly] structured text)]" << endl;
		cerr << "\t[--o.json (output results in JSON format)]" << endl;
		cerr << "\t[-d <max degeneracy> (default is " << DEFAULT_DEGEN << ")]" << endl;
		cerr << "\t[--trial <number of trials (default is " << DEFAULT_NUM_TRIAL << ")>]" << endl;
		cerr << "\t[--seed <random number seed> (default is time-based)]" << endl;
		cerr << "\t[--thread <maximum number of OpenMP threads> (default is all)]" << endl;
		cerr << "\t[--salt <salt concentration> (default is " << DEFAULT_SALT << ")" << endl;
		cerr << "\t[--primer.hairpin <max oligo hairpin Tm> (default is " << DEFAULT_MAX_HAIRPIN << ")" << endl;
		cerr << "\t[--primer.dimer <max oligo hetero/homo-dimer Tm> (default is " << DEFAULT_MAX_DIMER << ")" << endl;
		cerr << "\t[--count <total number of amplicons to produce> (default is " << DEFAULT_NUM_ASSAY << ")]" << endl;
		cerr << "\t[--optimize.top-down Search using maximally degenerate inital primer oligos (default is bottom up)]" << endl;
		//cerr << "\t[--multiplex (choose multiplex compatible assays)]" << endl;
		cerr << "\t[--optimize.5 (enable 5' primer search to optimize coverage)]" << endl;
		cerr << "\t[--no-optimize.5 (disable 5' primer search to optimize coverage; default)]" << endl;
		cerr << "\t[--optimize.3 (enable 3' primer search to optimize coverage)]" << endl;
		cerr << "\t[--no-optimize.3 (disable 3' primer search to optimize coverage; default)]" << endl;
		cerr << "\t[-v <verbosity level: silent, verbose, everything>]" << endl;
		cerr << "\t[--target.amplicon.min <minimum amplicon length> (default is " << DEFAULT_MIN_TARGET_AMPLICON << ")]" << endl;
		cerr << "\t[--target.amplicon.max <maximum amplicon length> (default is " << DEFAULT_MAX_TARGET_AMPLICON << ")]" << endl;
		cerr << "\t[--target.threshold <target detection threshold> (default is " << DEFAULT_TARGET_THRESHOLD << ")]" << endl;
		cerr << "\t[--target.search <target search multiplier> (default is " << DEFAULT_SEARCH_THRESHOLD_MULTIPLIER << ")]" << endl;
		cerr << "\t[--target.cover <minimum per primer pair coverage> (default is " << DEFAULT_MIN_TARGET_COVER << ")]" << endl;
		cerr << "\t[--target.ignore <defline key word to exclude a sequence>]" << endl;
		cerr << "\t[--target.normalize (normalize target weights per fasta file)]" << endl;
		cerr << "\t[--target.size.min <length> (minimum input target length in bp)]" << endl;
		cerr << "\t[--target.size.max <length> (maximum input target length in bp)]" << endl;
		cerr << "\t[--background.amplicon.min <minimum amplicon length>]" << endl;
		cerr << "\t[--background.amplicon.max <maximum amplicon length>]" << endl;
		cerr << "\t[--background.threshold <background detection threshold> (default is " << DEFAULT_BACKGROUND_THRESHOLD << ")]" << endl;
		cerr << "\t[--background.search <background search multiplier> (default is " << DEFAULT_SEARCH_THRESHOLD_MULTIPLIER << ")]" << endl;
		cerr << "\t[--background.cover <maximum per primer pair coverage> (default is " << DEFAULT_MAX_BACKGROUND_COVER << ")]" << endl;
		cerr << "\t[--background.ignore <defline key word to exclude a sequence>]" << endl;
		cerr << "\t[--background.normalize (normalize background weights per fasta file)]" << endl;
		cerr << "\t[--background.size.min <length> (minimum input background length in bp)]" << endl;
		cerr << "\t[--background.size.max <length> (maximum input background length in bp)]" << endl;
		cerr << "\t[--primer.size.min <minimum primer length> (default is " << DEFAULT_MIN_PRIMER << ")]" << endl;
		cerr << "\t[--primer.size.max <maximum primer length> (default is " << DEFAULT_MAX_PRIMER << ")]" << endl;
		cerr << "\t[--primer.tm.min <minimum primer melting temperature> (default is " << DEFAULT_MIN_PRIMER_TM << ")]" << endl;
		cerr << "\t[--primer.tm.max <maximum primer melting temperature> (default is " << DEFAULT_MAX_PRIMER_TM << ")]" << endl;
		cerr << "\t[--primer.strand <primer strand concentration> (default is " << DEFAULT_PRIMER_STRAND << ")]" << endl;
		cerr << "\t[--primer.taq-mama (use Taq MAMA rules for terminal primer mismatches; default is false)]" << endl;
		cerr << "\t[--pack.degen.max <max degen when packing words> (default is " << DEFAULT_PACK_MAX_DEGEN << ")]" << endl;
		cerr << "\t[--pack.gc.max <max fractional GC content when packing words> (default is " << DEFAULT_PACK_MAX_GC << ", disabled)]" << endl;
		cerr << "\t[--pack.gc.min <min fractional GC content when packing words> (default is " << DEFAULT_PACK_MIN_GC << ", disabled)]" << endl;
		cerr << "\t[--input.prefix <directory prefix for both target and background input genomes>]" << endl;

		// JSON-formatted input file is currently not recommended
		//cerr << "\t[--json <JSON file of configuration values>]" << endl;
		//cerr << "\t[--json.root <root key for JSON configuration file> (\"xxx|yyy|zzz|...\")]" << endl;
	
		quit = true;
		return;
	}

	if( !json_file.empty() ){
	
		if( !parse_json(json_file, json_root, 
			target_dir, target_dir_prefix,
			background_dir, background_dir_prefix) ){
			
			cerr << "Error parsing JSON file \"" << json_file << "\" with root key \"" 
				<< json_root << '"' << endl;
			quit = true;
			
			return;
		}
	}
	
	if( target_filename.empty() && target_dir.empty() ){

		cerr << "Please specify one or more target filenames (-t) or directories (-T)" << endl;
		
		quit = true;
		return;
	}

	if( output_filename.empty() ){

		cerr << "Please specify an output filename (-o)" << endl;
		
		quit = true;
		return;
	}
	
	if(degen == 0){
	
		cerr << "Please specify a valid maximum degeneracy (-d)" << endl;
		
		quit = true;
		return;
	}
	
	if(primer_range.second > WORD_LENGTH){
	
		cerr << "The maximum primer length must be <= " << WORD_LENGTH << endl;
		
		quit = true;
		return;
	}
		
	if(primer_range.second < primer_range.first){
	
		cerr << "The maximum primer length must be >= minimum primer length" << endl;
		
		quit = true;
		return;
	}
		
	if(primer_tm_range.first > primer_tm_range.second){
	
		cerr << "The maximum primer melting temperature must be >= minimum primer melting temperature" << endl;
		
		quit = true;
		return;
	}
		
	if(target_amplicon_range.second < target_amplicon_range.first){
	
		cerr << "The maximum target amplicon length must be >= minimum target amplicon length" << endl;
		
		quit = true;
		return;
	}
	
	if(background_amplicon_range.second < background_amplicon_range.first){
	
		cerr << "The maximum background amplicon length must be >= minimum background amplicon length" << endl;
		
		quit = true;
		return;
	}
	
	if(target_length_range.second < target_length_range.first){
	
		cerr << "The maximum target input sequence length must be >= minimum target input sequence length" << endl;
		
		quit = true;
		return;
	}
	
	if(background_length_range.second < background_length_range.first){
	
		cerr << "The maximum background input sequence length must be >= minimum background input sequence length" << endl;
		
		quit = true;
		return;
	}
	
	if(pack_max_gc < pack_min_gc){
	
		cerr << "The maximum packing GC content must be >= the minimum packing GC content" << endl;
		
		quit = true;
		return;
	}
	
	if(seed == 0){
		seed = time(NULL);
	}
	
	if(salt < 0.0f){
	
		cerr << "Please specify a salt concentration > 0" << endl;
		
		quit = true;
		return;
	}
	
	if(primer_strand < 0.0f){
	
		cerr << "Please specify a [primer strand] concentration > 0" << endl;
		
		quit = true;
		return;
	}
		
	// Make sure that we don't attempt to read the same data more than once
	make_set(target_filename);
	make_set(background_filename);
	
	// Strip any trailing path separators from the directory names to ensure
	// that group naming is consistent
	strip_trailing_path_separtor(target_dir);
	strip_trailing_path_separtor(background_dir);
	
	make_set(target_dir);
	make_set(background_dir);
	
	// If the user has set the directory prefix, but not the more specific
	// target or background directory prefixs, set them here.
	if( target_dir_prefix.empty() ){
		target_dir_prefix = dir_prefix;
	}
	
	if( background_dir_prefix.empty() ){
		background_dir_prefix = dir_prefix;
	}
	
	// Recursively search each input directory or file and extract all subdirectories that contain
	// fasta files. Each subdirectory (or file name) becomes a group name for the fasta files contained
	// in the subdirectory (or file). All of the fasta records associated with a group will be
	// concatinated into a single sequence.
	for(deque<string>::const_iterator i = target_dir.begin();i != target_dir.end();++i){
	
		string target_path;
		
		if( target_dir_prefix.empty() ){
			target_path = *i;
		}
		else{
			target_path = target_dir_prefix + PATH_SEPARATOR + *i;
		}
		
		if(find_groups(target_path, target_groups) == false){
		
			cerr << "Invalid target path: " << target_path << endl;
			
			quit = true;
			return;
		}
	}
	
	for(deque<string>::const_iterator i = background_dir.begin();i != background_dir.end();++i){
		
		string background_path;
		
		if( background_dir_prefix.empty() ){
			background_path = *i;
		}
		else{
			background_path = background_dir_prefix + PATH_SEPARATOR + *i;
		}
		
		if(find_groups(background_path, background_groups) == false){
			
			cerr << "Invalid background path: " << *i << endl;
			
			quit = true;
			return;
		}
	}
};


Options::Verbosity parse_verbosity(string m_buffer)
{
	// Make buffer lower case
	for(string::iterator i = m_buffer.begin();i != m_buffer.end();++i){
		*i = tolower(*i);
	}
	
	if(m_buffer == "silent"){
		return Options::SILENT;
	}
	
	if(m_buffer == "verbose"){
		return Options::VERBOSE;
	}
	
	if(m_buffer == "everything"){
		return Options::EVERYTHING;
	}
	
	return Options::UNKNOWN_VERBOSITY;
}

bool Options::parse_json(const string &m_filename, const string &m_root_key, 
	deque<string> &m_target_dir, string &m_target_dir_prefix,
	deque<string> &m_background_dir, string &m_background_dir_prefix)
{
	try{
		ifstream fin( m_filename.c_str() );

		if(!fin){
			throw __FILE__ ":Options::parse_json: Unable to open JSON file for reading";
		}
		
		string line;
		string buffer;
		
		while( getline(fin, line) ){
			buffer += line;
		}
		
		fin.close();
		
		json::JSON conf(buffer);
		
		const deque<string> root = parse_keys(m_root_key);
		
		if( !root.empty() ){
			
			// Change the root
			for(deque<string>::const_iterator i = root.begin();i != root.end();++i){
				
				if( (conf.get_type() != json::JSON::JSON_MAP) || !conf.has_key(*i) ){
					
					cerr << "Could not find JSON key: \"" << *i << '"' << endl;
					throw __FILE__ ":Options::parse_json: Unable to find JSON key";
				}
				
				conf = conf[*i];
			}
		}
		
		if(conf.get_type() != json::JSON::JSON_MAP){
			throw __FILE__ ":Options::parse_json: Specified root key does not yeild map";
		}
		
		/////////////////////////////////////////////////////////////////////////////////////
		// Look for valid configuration key/value pairs
		/////////////////////////////////////////////////////////////////////////////////////
		
		if( conf.has_key("output_file") ){
			output_filename = conf["output_file"].get_string();
		}
		
		if( conf.has_key("target_species") ){
			
			// This key needs to point to an array
			if(conf["target_species"].get_type() != json::JSON::JSON_ARRAY){
				throw __FILE__ ":Options::parse_json: \"target_species\" key does not yeild array";
			}
			
			const size_t num_target = conf["target_species"].size();
			
			for(size_t i = 0;i < num_target;++i){
				
				if( !conf["target_species"][i].has_key("value") ){
					
					cerr << "Unable to find conf[\"target_species\"][" << i << "][\"value\"]" << endl;
					throw __FILE__ ":Options::parse_json: Missing \"value\" key";
				}
				
				m_target_dir.push_back( conf["target_species"][i]["value"].get_string() );
				
				// For now, only replace white space with '_'
				m_target_dir.back() = replace_special_with( '_', m_target_dir.back() );
			}
		}
		
		if( conf.has_key("bg_species") ){
			
			// This key needs to point to an array
			if(conf["bg_species"].get_type() != json::JSON::JSON_ARRAY){
				throw __FILE__ ":Options::parse_json: \"bg_species\" key does not yeild array";
			}
			
			const size_t num_background = conf["bg_species"].size();
						
			for(size_t i = 0;i < num_background;++i){
				
				if( !conf["bg_species"][i].has_key("value") ){
					
					cerr << "Unable to find conf[\"bg_species\"][" << i << "][\"value\"]" << endl;
					throw __FILE__ ":Options::parse_json: Missing \"value\" key";
				}
				
				m_background_dir.push_back( conf["bg_species"][i]["value"].get_string() );
				
				// For now, only replace white space with '_'
				m_background_dir.back() = replace_special_with( '_', m_background_dir.back() );
			}
		}
		
		if( conf.has_key("maxDeg") ){			
			degen = conf["maxDeg"].get_number(true);
		}
		
		if( conf.has_key("numTrial") ){			
			num_trial = conf["numTrial"].get_number(true);
		}
		
		// Multiplex assay design is the default
		//if( conf.has_key("multiplex") ){
		//	use_multiplex = conf["multiplex"].get_bool(true);
		//}
		
		if( conf.has_key("target_min_amplicon") ){
			target_amplicon_range.first = conf["target_min_amplicon"].get_number(true);
		}
		
		if( conf.has_key("target_max_amplicon") ){
			target_amplicon_range.second = conf["target_max_amplicon"].get_number(true);
		}
		
		if( conf.has_key("target_detect_threshold") ){
			target_threshold = conf["target_detect_threshold"].get_number(true);
		}
		
		if( conf.has_key("bg_min_amplicon") ){
			background_amplicon_range.first = conf["bg_min_amplicon"].get_number(true);
		}
		
		if( conf.has_key("bg_max_amplicon") ){
			background_amplicon_range.second = conf["bg_max_amplicon"].get_number(true);
		}
		
		if( conf.has_key("bg_detect_threshold") ){
			background_threshold = conf["bg_detect_threshold"].get_number(true);
		}
		
		if( conf.has_key("min_primer_len") ){
			primer_range.first = conf["min_primer_len"].get_number(true);
		}
		
		if( conf.has_key("max_primer_len") ){
			primer_range.second = conf["max_primer_len"].get_number(true);
		}
		
		if( conf.has_key("min_primer_tm") ){
			primer_tm_range.first = conf["min_primer_tm"].get_number(true);
		}
		
		if( conf.has_key("max_primer_tm") ){
			primer_tm_range.second = conf["max_primer_tm"].get_number(true);
		}
		
		if( conf.has_key("primer_strandcon") ){
			primer_strand = conf["primer_strandcon"].get_number(true);
		}
		
		if( conf.has_key("salt") ){
			salt = conf["salt"].get_number(true);
		}
		
		if( conf.has_key("hairpin_tm") ){
			max_hairpin = conf["hairpin_tm"].get_number(true);
		}
		
		if( conf.has_key("dimer_tm") ){
			max_dimer = conf["dimer_tm"].get_number(true);
		}
		
		if( conf.has_key("count") ){
			num_assay = conf["count"].get_number(true);
		}
		
		if( conf.has_key("optimize5") ){
			optimize_5 = conf["optimize5"].get_bool(true);
		}
		
		if( conf.has_key("optimize3") ){
			optimize_3 = conf["optimize3"].get_bool(true);
		}
		
		if( conf.has_key("target_search_factor") ){
			target_search_multiplier = conf["target_search_factor"].get_number(true);
		}
		
		if( conf.has_key("target_min_cov") ){
			min_target_cover = conf["target_min_cov"].get_number(true);
		}
		
		if( conf.has_key("target_normalize") ){
			normalize_target_weight_per_file = conf["target_normalize"].get_bool(true);
		}
		
		if( conf.has_key("target_minLen") ){
			target_length_range.first = conf["target_minLen"].get_number(true);
		}
		
		if( conf.has_key("target_maxLen") ){
			target_length_range.second = conf["target_maxLen"].get_number(true);
		}
		
		if( conf.has_key("bg_search_factor") ){
			background_search_multiplier = conf["bg_search_factor"].get_number(true);
		}
		
		if( conf.has_key("bg_max_cov") ){
			max_background_cover = conf["bg_max_cov"].get_number(true);
		}
		
		if( conf.has_key("bg_normalize") ){
			normalize_background_weight_per_file = conf["bg_normalize"].get_bool(true);
		}
		
		if( conf.has_key("bg_minLen") ){
			background_length_range.first = conf["bg_minLen"].get_number(true);
		}
		
		if( conf.has_key("bg_maxLen") ){
			background_length_range.second = conf["bg_maxLen"].get_number(true);
		}
		
		if( conf.has_key("primer_taq_mama") ){
			use_taq_mama = conf["primer_taq_mama"].get_bool(true);
		}
				
		if( conf.has_key("seed") ){
			seed = conf["seed"].get_number(true);
		}
		
		if( conf.has_key("max_pack_degen") ){
			pack_max_degen = conf["max_pack_degen"].get_number(true);
		}
		
		if( conf.has_key("min_pack_gc") ){
			pack_min_gc = conf["min_pack_gc"].get_number(true);
		}
		
		if( conf.has_key("max_pack_gc") ){
			pack_max_gc = conf["max_pack_gc"].get_number(true);
		}
		
		if( conf.has_key("input_prefix") ){
			m_target_dir_prefix = m_background_dir_prefix = conf["input_prefix"].get_string();
		}
		
		if( conf.has_key("target_prefix") ){
			m_target_dir_prefix = conf["target_prefix"].get_string();
		}
		
		if( conf.has_key("background_prefix") ){
			m_background_dir_prefix = conf["background_prefix"].get_string();
		}
	}
	catch(const char *error){
	
		cerr << "Error reading JSON file: " << error << endl;
		return false;
	}
	catch(...){
	
		cerr << "Unknown error reading JSON file" << endl;
		return false;
	}
	
	return true;
}

string replace_special_with(char m_replacement, const string &m_str)
{
	string ret(m_str);
	
	for(string::iterator i = ret.begin();i != ret.end();++i){
		
		switch(*i){
			case ' ':
			case '\t':
			case '\n':
			case '\r':
				*i = m_replacement;
				break;
		};
	}
	
	return ret;
}

deque<string> parse_keys(const string &m_key_str)
{
	deque<string> ret;
	
	const char delim = '|';
	
	if( m_key_str.empty() ){
		return ret;
	}
	
	ret.push_back( string() );
	
	for(string::const_iterator i = m_key_str.begin();i != m_key_str.end();++i){
		
		if(*i == delim){
			ret.push_back( string() );
		}
		else{
			ret.back().push_back(*i);
		}
	}
	
	return ret;
}

string tolower(const string &m_str)
{
	string ret(m_str);
	
	for(string::iterator i = ret.begin();i != ret.end();++i){
		*i = tolower(*i);
	}
	
	return ret;
}

bool find_groups(const string &m_path, unordered_multimap<string, string> &m_group)
{	
	// Is the input path a directory?
	struct stat path_info;

	if(stat(m_path.c_str(), &path_info) != 0){
		return false;
	}
		
	// A file or directory with this name already exists
	if( S_ISDIR(path_info.st_mode) ){ // A directory
		
		// Read the contents of the directory
		DIR *dp;
		struct dirent *dir;
		
		dp = opendir( m_path.c_str() );

		if(dp == NULL){
			return false;
		}
		
		deque<string> dir_to_read;
		
		while( ( dir = readdir(dp) ) ){
		
			// Skip any removed files or diretories
			if(dir->d_ino == 0){
				continue;
			}

			// Skip the special directories "." and ".."
			if(strcmp(dir->d_name, ".") == 0){
				continue;
			}

			if(strcmp(dir->d_name, "..") == 0){
				continue;
			}
			
			const string name = m_path + PATH_SEPARATOR + dir->d_name;
				
			// Is this entry a directory or file name?
			struct stat dir_info;

			stat(name.c_str(), &dir_info);

			if( S_ISDIR(dir_info.st_mode) ){
			
				// Save the names of the subdirectories to search. Don't
				// search them now to avoid the risk of opening too many
				// directories for reading (since we do not know the system
				// limit).
				dir_to_read.push_back(name);
			}
			else if( S_ISREG(dir_info.st_mode) && 
				find_file_extension(name, allowed_fasta_extions) ) {
				
				m_group.insert( make_pair(m_path, name) );
			}
		}

		closedir(dp);
		
		for(deque<string>::const_iterator i = dir_to_read.begin();i != dir_to_read.end();++i){
			
			if(find_groups(*i, m_group) == false){
				return false;
			}
		}
		
		return true;
		
	}
	else if( S_ISREG(path_info.st_mode) ) { // A regular file
		
		// Is this a fasta file?
		if( find_file_extension(m_path, allowed_fasta_extions) ){
		
			m_group.insert( make_pair(m_path, m_path) );
			
			return true;
		}
		
		// It is an error of the user has specified an file that does
		// not have a valid fasta extension
		return false;
	}
	
	// The user has provide a path that does not point to a valid directory of file
	return false;
}

bool find_file_extension(const string &m_path, const char** m_ext)
{
	
	for(const char** ptr = m_ext;*ptr != NULL;++ptr){
		
		size_t loc = m_path.find(*ptr);
		
		if(loc != string::npos){
			
			if( (strlen(*ptr) + loc) == m_path.size() ){
				return true;
			}
		}
	}
	
	return false;
}

void strip_trailing_path_separtor(deque<string> &m_str)
{
	for(deque<string>::iterator i = m_str.begin();i != m_str.end();++i){
		
		// Some older C++ compilers do not have string::pop_back()
		const size_t init_len = i->size();
		size_t len = init_len;
		
		while( (len > 0) && ( (*i)[len - 1] == PATH_SEPARATOR ) ){
			--len;
		}
		
		if(len != init_len){
			*i = i->substr(0, len);
		}
	}
}
