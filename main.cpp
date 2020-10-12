// PCRamp: PCR-based target enrichment for sequencing
// J. D. Gans
// Bioscience Division, B-10
// Los Alamos National Laboratory
// Thursday Feb 7, 2020

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <deque>
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

// Still no openmp support on the clang compiler that ships with OSX
#ifdef _OPENMP
#include <omp.h>
#endif // _OPENM

#include "update.h"
#include "pcramp.h"
#include "assay.h"

using namespace std;

// Global variables for MPI
int mpi_numtasks;
int mpi_rank;

// MPI messages
#define	BEST_ASSAY	1000
#define	BEST_ASSAY_DATA	1001

string time_to_str(const time_t &m_elapsed);
void sequence_summary(const string &m_prefix, const deque<Sequence> &m_seq, ostream &m_out,
	const Options::OutputFormat &m_format);
float weighted_coverage(const BitSet &m_match, const deque<Sequence> &m_seq);
void reduce_best_assay(PCR &m_assay, Score &m_score, BitSet &m_target_match, BitSet &m_background_match,
	deque<Sequence> &m_amplicons, deque<AmpliconBounds> &m_bounds);
string truncate_prefix(const string &m_str, const size_t &m_max_len);

int main(int argc, char *argv[])
{
	try{
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

		time_t profile = time(NULL);

		const bool is_root = (mpi_rank == 0);
		
		Options opt;
				
		if(is_root){
		
			opt.load(argc, argv);
			
			// Partition the number of trials across all MPI tasks
			opt.num_trial = max(1u, opt.num_trial/mpi_numtasks + ( (opt.num_trial%mpi_numtasks == 0) ? 0u : 1u) );
		}
		
		// Share the options with all workers
		broadcast(opt, mpi_rank, 0);
		
		if(opt.quit){
			
			MPI_Finalize();
			return EXIT_SUCCESS;
		}
		
		// Define the allowed optimization and relaxation moves
		vector<Move> optimization_moves;
		
		if(opt.degen > 1){
		
			optimization_moves.push_back(IncreaseDegeneracy);
			optimization_moves.push_back(DecreaseDegeneracy);
		}
		
		if(opt.optimize_5){
			
			optimization_moves.push_back(Trim5);
			optimization_moves.push_back(Grow5);
		}
		
		if(opt.optimize_3){
			
			optimization_moves.push_back(Trim3);
			optimization_moves.push_back(Grow3);
		}
		
		// Make sure that all of the random number generators get initialized.		
		MPI_Bcast( (void*)&(opt.seed), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
		
		// Every worker gets a unique seed by adding its rank to the global seed (which is
		// used by the rank 0 task).
		opt.seed += mpi_rank;
		
		if(opt.max_thread > 0){
		
			#ifdef _OPENMP
			omp_set_num_threads(opt.max_thread);
			#endif // _OPENMP
		}
		
		unsigned int global_seed = opt.seed;
		
		srand(global_seed);
		srandom(global_seed);
		
		ofstream fnull("/dev/null");
		
		if(!fnull){
			throw "Unable to open /dev/null";
		}
		
		ostream &vout = ( (opt.output_filter == Options::SILENT) || !is_root ) ? fnull : cerr;
		
		ofstream fout(is_root ? opt.output_filename.c_str() : "/dev/null");
		
		if(!fout){
			throw __FILE__ ": Unable to open output file for writing";
		}
		
		switch(opt.output_format){
			case Options::TEXT_OUTPUT:
				fout << "PCRamp version " 
					<< PCRAMP_MAJOR_VERSION << '.' 
					<< PCRAMP_MINOR_VERSION << endl;

				// Write the command line arguments to disk
				fout << "Command line:";

				for(int i = 0;i < argc;++i){
					fout << ' ' << argv[i];
				}

				fout << endl;

				fout << "Random number seed = " << opt.seed << endl;
				break;
			case Options::JSON_OUTPUT:
				fout << "{\n\t\"program\":\"PCRamp\",\n"
					<< "\t\"version\":\"" << PCRAMP_MAJOR_VERSION << '.'
					<< PCRAMP_MINOR_VERSION << "\",\n\t"
					<< "\"command line\":\"" << argv[0];
				
				for(int i = 1;i < argc;++i){
					fout << ' ' << argv[i];
				}

				fout << "\",\n\t\"seed\":" << opt.seed << ',' << endl;
				
				break;
			default:
				throw __FILE__ ":main: Unknown output format";
		};
		
		vout << "PCRamp version " 
			<< PCRAMP_MAJOR_VERSION << "." 
			<< PCRAMP_MINOR_VERSION << endl;
		
		if(mpi_numtasks > 1){
			vout << "Running in MPI mode with " << mpi_numtasks << " workers" << endl;
		}
		
		vout << "Random number seed = " << global_seed << endl;
		vout << "background.threshold = " << opt.background_threshold << endl;
		vout << "target.threshold = " << opt.target_threshold << endl;
		vout << "Using a " << (opt.top_down_search ? "top-down" : "bottom-up") 
			<< " search strategy" << endl;
		vout << "[Na+] = " << opt.salt << endl;
		vout << "Max oligo hairpin Tm <= " << opt.max_hairpin << endl;
		vout << "Allowed input target sequence length range = [" 
			<< opt.target_length_range.first << ", "
			<< opt.target_length_range.second << "]" << endl;
		vout << "Allowed input background sequence length range = [" 
			<< opt.background_length_range.first << ", "
			<< opt.background_length_range.second << "]" << endl;
		
		vout << "Local oligo optimization moves:" << endl;
		
		if(opt.degen > 1){
			vout << "\tDegeneracy allowed" << endl;
		}
		else{
			vout << "\tDegeneracy not allowed" << endl;
		}
		
		if(opt.optimize_5){
			vout << "\t5' oligo search allowed" << endl;
		}
		else{
			vout << "\t5' oligo search not allowed" << endl;
		}
		
		if(opt.optimize_3){
			vout << "\t3' oligo search allowed" << endl;
		}
		else{
			vout << "\t3' oligo search not allowed" << endl;
		}
		
		// When reading multi-fasta record targets or backgrounds, how many "virtual" bases will
		// separated distinct sequences?
		const size_t target_group_padding = 1;
		const size_t background_group_padding = 1;
		
		vout << '\t' << opt.target_amplicon_range.first << " <= target amplicon length <= " 
			<< opt.target_amplicon_range.second << endl;
		vout << '\t' << opt.background_amplicon_range.first << " <= background amplicon length <= " 
			<< opt.background_amplicon_range.second << endl;
		vout << '\t' << opt.primer_range.first << " <= primer length <= " 
			<< opt.primer_range.second << endl;
		vout << '\t' << opt.primer_tm_range.first << " <= primer Tm <= " 
			<< opt.primer_tm_range.second << endl;
			
		if(opt.use_taq_mama){
			vout << "** Using Taq-MAMA primer binding rules **" << endl;
		}
		
		if(opt.use_multiplex){
			vout << "** Designing multiplex-compatible assays **" << endl;
		}
			
		vout << "Maximum oligo degeneracy = " << opt.degen << endl;
		
		if(mpi_numtasks > 1){
		
			vout << "Number of assay design trials per MPI task = " << opt.num_trial << endl;
			vout << "Total number of assay design trials = " << opt.num_trial*mpi_numtasks << endl;
		}
		else{
			vout << "Number of assay design trials = " << opt.num_trial << endl;
		}
		
		// Read the user-supplied sequences
		deque<Sequence> target_seq;
		deque<Sequence> background_seq;
		
		// The multiplex background sequences are only used for multipex PCR design.
		// These sequences are *not* user specified, but store all of the current *amplicon* (minus primer binding
		// site) sequences. These amplicon sequences are incermentally updated during the greedy assay design process
		// an are used to insure that *new* primer designs do not overlap exising amplicons (which could create
		// spurious amplicons with the potential to inhibit intended PCR amplifications).
		deque<Sequence> multiplex_background_seq;
		
		if(is_root){
		
			// Each fasta record in opt.target_filename is a distinct detection target
			for(deque<string>::const_iterator i = opt.target_filename.begin();i != opt.target_filename.end();++i){
				
				const size_t seq_begin = target_seq.size();
				
				parse_fasta(*i, target_seq, 
					max(opt.target_amplicon_range.first, opt.target_length_range.first),
					opt.target_length_range.second,
					opt.target_ignore);
				
				const size_t seq_end = target_seq.size();
				
				if(opt.normalize_target_weight_per_file){

					const double w = (seq_end == seq_begin) ? 1.0 : 1.0/(seq_end - seq_begin);

					for(size_t j = seq_begin;j < seq_end;++j){
						target_seq[j].weight(w);
					}
				}
			}

			const vector<string> groups = keys(opt.target_groups);
			
			if( !groups.empty() ){
				
				UpdateInfo info("Reading target groups: ", vout);
				
				for(vector<string>::const_iterator i = groups.begin();i != groups.end();++i){

					if( ignore_record(*i, opt.target_ignore) ){
						
						info << "skipping " << truncate_prefix(*i, 50);
						info.flush();
					
						continue;
					}
					
					info << ( 100.0*( 1.0 + ( i - groups.begin() ) ) )/groups.size() 
						<< "% " << truncate_prefix(*i, 50);
					info.flush();
						
					target_seq.push_back( Sequence() );

					Sequence &ref = target_seq.back();

					// Use the group name as the define ...
					string name = *i;

					// ... but strip off the target dir prefix (if it
					// is defined). Also strip off any leading '/' left in the name
					if(name.find(opt.target_dir_prefix) == 0){

						name = name.substr( opt.target_dir_prefix.size(), 
							name.size() - opt.target_dir_prefix.size() );
						
						while(!name.empty() && (name[0] == PATH_SEPARATOR) ){
							name = name.substr(1, name.size() - 1);
						}
					}

					ref.defline(name); 

					ref.active(true);

					typedef unordered_multimap<string, string>::const_iterator I;

					const pair<I, I> range = opt.target_groups.equal_range(*i);

					for(I j = range.first;j != range.second;++j){

						append_fasta_group(j->second, ref,
							max(opt.target_amplicon_range.first, opt.target_length_range.first), 
							opt.target_length_range.second,
							target_group_padding, opt.target_ignore);
					}
					
					// Don't save zero length sequences that were skipped due to the sequence length
					// constraints
					if( ref.empty() ){
						target_seq.pop_back();
					}
				}
				
				info.close();
			}
		}

		broadcast(target_seq, mpi_rank, 0);
		
		if(is_root){
			
			for(deque<string>::const_iterator i = opt.background_filename.begin();i != opt.background_filename.end();++i){
				
				const size_t seq_begin = background_seq.size();
				
				parse_fasta(*i, background_seq, 
					max(opt.background_amplicon_range.first, opt.background_length_range.first),
					opt.background_length_range.second,
					opt.background_ignore);
				
				const size_t seq_end = background_seq.size();
				
				if(opt.normalize_background_weight_per_file){

					const double w = (seq_end == seq_begin) ? 1.0 : 1.0/(seq_end - seq_begin);

					for(size_t j = seq_begin;j < seq_end;++j){
						background_seq[j].weight(w);
					}
				}
			}
			
			const vector<string> groups = keys(opt.background_groups);
			
			if( !groups.empty() ){

				UpdateInfo info("Reading background groups: ", vout);
				
				for(vector<string>::const_iterator i = groups.begin();i != groups.end();++i){

					if( ignore_record(*i, opt.background_ignore) ){
						
						info << "skipping " << truncate_prefix(*i, 50);
						info.flush();
						
						continue;
					}

					info << ( 100.0*( 1.0 + ( i - groups.begin() ) ) )/groups.size() 
						<< "% " << truncate_prefix(*i, 50);
					info.flush();
					
					background_seq.push_back( Sequence() );

					Sequence &ref = background_seq.back();

					// Use the group name as the define ...
					string name = *i;

					// ... but strip off the background dir prefix (if it
					// is defined). Also strip off any leading '/' left in the name
					if(name.find(opt.background_dir_prefix) == 0){

						name = name.substr( opt.background_dir_prefix.size(), 
							name.size() - opt.background_dir_prefix.size() );
						
						while(!name.empty() && (name[0] == PATH_SEPARATOR) ){
							name = name.substr(1, name.size() - 1);
						}
					}

					ref.defline(name);
					ref.active(true);

					typedef unordered_multimap<string, string>::const_iterator I;

					const pair<I, I> range = opt.background_groups.equal_range(*i);

					for(I j = range.first;j != range.second;++j){

						append_fasta_group(j->second, ref,
							max(opt.background_amplicon_range.first, opt.background_length_range.first), 
							opt.background_length_range.second,
							background_group_padding, opt.background_ignore);
					}
					
					// Don't save zero length sequences that were skipped due to the sequence length
					// constraints
					if( ref.empty() ){
						background_seq.pop_back();
					}
				}
				
				info.close();
			}
		}
		
		broadcast(background_seq, mpi_rank, 0);
		
		const unsigned int num_target_seq = target_seq.size();
		const unsigned int num_background_seq = background_seq.size();
			
		// Report sequence statistics
		sequence_summary("target sequence summary", target_seq, fout, opt.output_format);
		sequence_summary("Target:", target_seq, vout, Options::TEXT_OUTPUT);
		
		sequence_summary("background sequence summary", background_seq, fout, opt.output_format);
		sequence_summary("Background:", background_seq, vout, Options::TEXT_OUTPUT);
		
		// Store the best assays found during the entire search
		deque<PCR> assay_pool;
		deque<BitSet> pool_background; // Store the background coverage of each assay in the pool
		
		// The multiplex backgrounds are only used for multiplex PCR-bases assay design, and are incrementally
		// updated to store the existing amplicon sequences (to insure that new primers do *not* bind to existing
		// amplicons).
		MULTIMAP<Word, WordMatch> multiplex_background_db;
		vector<Word> multiplex_background_keys;

		/////////////////////////////////////////////////////////////////////////////////////////////////////		

		unsigned int assay_iteration = 0;

		if(opt.output_format == Options::JSON_OUTPUT){
			fout << "\t\"assays\":[\n";
		}
		
		// Once all targets have been detected, we increment the major assay id
		// an reset the minor assay id. This process continues until we have
		// designed the requested number of assays or we can not obtain an
		// additional assay design.
		unsigned int major_assay_id = 1;
		unsigned int minor_assay_id = 1;

		while(true){ // Keep looping until we have designed num_assay assays

			++assay_iteration;
			
			// Count the number of remaining targets
			unsigned int targets_remaining = 0;

			for(unsigned int i = 0;i < num_target_seq;++i){

				if( target_seq[i].active() ){
					++targets_remaining;
				}
			}

			if(opt.output_filter > Options::VERBOSE){
				vout << '\t' << targets_remaining << " targets remaining in current iterations" << endl;
			}

			if(targets_remaining == 0){
				
				// We have detected all targets! Reset the active flag for
				// all targets so that we can design the requested number
				// of assays
				for(unsigned int i = 0;i < num_target_seq;++i){
					target_seq[i].active(true);
				}

				targets_remaining = num_target_seq;

				++major_assay_id;
				minor_assay_id = 1;
			}
			
			if(opt.output_format == Options::JSON_OUTPUT){
				
				if(assay_iteration > 1){
					fout << ",\n";
				}
				
				fout << "\t\t{\n\t\t\t\"id\":" << major_assay_id << '.' << minor_assay_id << ",\n";
			}
			
			vout << "Design iteration " << assay_iteration << endl;
			
			if(opt.output_format == Options::TEXT_OUTPUT){
			
				fout << "###########################################################################################" << endl;
				fout << "# Attempting to detect " << targets_remaining << " remaining targets" << endl;
			}
			
			PCR best_assay;
			BitSet best_background_match(num_background_seq, false);
			Score best_score;
			
			vector<PCR> trial_assays(opt.num_trial);
			
			// Randomly sample all of the trial assays *before* we attempt to optimize them
			#pragma omp parallel
			{

				// A melting temperature engine that will be used to filter randomly generated assays
				// and select valid assays
				NucCruc melt;

				melt.salt(opt.salt);

				// Generate a unique local seed for each thread
				unsigned int local_seed = 0;

				#pragma omp critical
				local_seed = rand_r(&global_seed);

				#pragma omp for
				for(unsigned int t = 0;t < opt.num_trial;++t){

					// Randomly select a target and generate an assay
					trial_assays[t].random_assay(target_seq, 
						melt, opt, local_seed, vout);
						
					
					#ifdef TARGET_SPIKE_IN
					trial_assays[t].oligo(FORWARD, Word("AGAAGGCTCGCCAAAATAAACG") );
					trial_assays[t].oligo(REVERSE, Word("TTGGACACACAAAAAAGAA") );
					
					trial_assays[t].center();
					#endif // TARGET_SPIKE_IN
				}
			}
			
			// Only pack the parts of the background sequences that are potential matches for one or
			// more trial assays.
			MULTIMAP<Word, WordMatch> background_db;
			vector<Word> background_keys;

			size_t num_active_background = 0;
			float active_background_norm = 0.0f;
			
			for(deque<Sequence>::const_iterator i = background_seq.begin();i != background_seq.end();++i){
				num_active_background += ( i->active() ? 1 : 0 );
			}
				
			time_t profile = time(NULL);
			
			if(num_background_seq > 0){

				UpdateInfo info("\tPreparing background for search: ", vout);
				const unsigned int update_every = max(1.0f, num_active_background*0.01f);

				for(unsigned int i = 0;i < num_background_seq;++i){

					if( !background_seq[i].active() ){
						continue;
					}
					
					MULTIMAP<Word, WordMatch> local_db;
					
					// Include slightly shorter background words (as low as the 90% of the
					// specified min_oligo_length) to prevent assay oligos that target the
					// ends of sequences from falsely appearing to be unique due to the
					// background failing to contain the partially matching oligo with the
					// *correct* offset.
					background_seq[i].pack( local_db, i, 
						opt.pack_max_degen, 
						0.0, 1.0, // Don't G+C filter the background sequences
						opt.min_oligo_length()*0.9 );
					
					// Only include the parts of this background (currently in local_db) that
					// are potential matches to one of the candidate assays
					select_words(background_db, local_db, trial_assays,
						opt.optimize_5, opt.optimize_3,
						opt.background_threshold*opt.background_search_multiplier);
					
					active_background_norm += background_seq[i].weight();

					if( ( (i + 1) % update_every) == 0 ){

						info << (i + 1)*100.0f/num_background_seq << '%';
						info.flush();
					}
					
					// DEBUG
					//cerr << "[" << mpi_rank << "] |background_db| = " 
					//	<< background_db.size() << "; |local| = " 
					//	<< local_db.size() << endl;
				}

				// Make the background multimap iterable
				background_db.sort();
	
				info.close();

				profile = time(NULL) - profile;
			
				vout << "\t\tIndexed background in " << profile << " sec" << endl;

				background_keys = keys(background_db);

				vout << "\tBackground word table has " << background_db.size() << " entries" << endl;
				vout << "\tFound " << background_keys.size() << " unique word keys (" << ( 100.0f*background_keys.size()/background_db.size() ) 
					<< "% of entries)" << endl;
			}

			unsigned int num_active_target = 0;
			float active_target_norm = 0.0f;

			UpdateInfo info("\tPreparing targets for search: ", vout);
			unsigned int update_every = 
				max(1.0f, num_target_seq*0.01f);

			MULTIMAP<Word, WordMatch> target_db;

			profile = time(NULL);

			for(unsigned int i = 0;i < num_target_seq;++i){

				if( target_seq[i].active() ){

					++num_active_target;
					active_target_norm += target_seq[i].weight();
				}
				else{
					// Don't pack target sequences that
					// we have already covered
					continue;
				}

				MULTIMAP<Word, WordMatch> local_db;
				
				target_seq[i].pack( local_db, i, 
					opt.pack_max_degen, 
					opt.pack_min_gc, 
					opt.pack_max_gc,
					opt.min_oligo_length() );

				// Only include the parts of this target (currently in local_db) that
				// are potential matches to one of the candidate assays
				select_words(target_db, local_db, trial_assays,
					opt.optimize_5, opt.optimize_3,
					opt.target_threshold*opt.target_search_multiplier);
				
				if( ( (i + 1) % update_every) == 0 ){

					info << (i + 1)*100.0f/num_target_seq << '%';
					info.flush();
				}
			}
			
			// Make the background multimap iterable
			target_db.sort();
				
			info.close();
			
			profile = time(NULL) - profile;
			
			vout << "\t\tIndexed targets in " << profile << " sec" << endl;
			
			vout << "\t\tNumber of active target sequences = " << num_active_target 
				<< " (total weight = " << active_target_norm << ")" << endl;

			// Use a vector instead of a set for efficient access later ...
			const vector<Word> target_keys = keys(target_db);

			vout << "\tTarget word table has " << target_db.size() << " entries" << endl;
			vout << "\tFound " << target_keys.size() << " unique word keys (" << ( 100.0f*target_keys.size()/target_db.size() ) 
				<< "% of entries)" << endl;

			#pragma omp parallel
			{
				
				// A melting temperature engine that will be used to filter randomly generated assays
				// and select valid assays
				NucCruc melt;

				melt.salt(opt.salt);
				
				#pragma omp for
				for(unsigned int t = 0;t < opt.num_trial;++t){
					
					if(opt.top_down_search){

						// Build the maximimally degenerate assay by comparing
						// the initial random assay to similar target sequences and
						// taking the set union of all bases in all oligos until
						// we reach the maximim allowed degeneracy
						
						const bool valid = make_degenerate(trial_assays[t], target_keys, target_db, 
							target_seq, melt, opt, vout);

						if(!valid){
							continue;
						}

					} // bottom-up otherwise
			
					Score s = optimize(trial_assays[t], optimization_moves, 
						target_keys, target_db, target_seq, 
						background_keys, background_db, background_seq, 
						multiplex_background_keys, multiplex_background_db, multiplex_background_seq, 
						assay_pool, opt, vout);
					
					// If the approximate background coverage is too high, then this assay
					// is already invalid (since the approximate background coverage is a
					// lower bound to the true background coverage)
					if( (s.background_coverage > opt.max_background_cover) ||
					    (s.target_coverage < opt.min_target_cover) ){

						continue;
					}
										
					// Recompute the background coverage.
					s.background_coverage = 0;

					if(opt.use_multiplex){

						// Test the optimized assay for multiplex compatibility
						bool is_compatible = true;

						for(deque<PCR>::const_iterator i = assay_pool.begin();
							( i != assay_pool.end() ) && is_compatible;++i){

							is_compatible = i->multiplex_compatible(melt, opt, trial_assays[t]);
						}

						if(!is_compatible){

							continue;
						}

						if(best_score < s){

							// Find Smith-Waterman alignment-based amplicon matches
							BitSet multiplex_background_match(multiplex_background_seq.size(), false);

							// Only match the background sequences. Targets are NOT allowed to have gaps, so our gap free word
							// matching is good enough for counting target matches. The actual identity of target matches will
							// be computed as the end.
							trial_assays[t].find_multiplex_background_match(multiplex_background_match, 
								multiplex_background_seq, opt, vout);

							s.background_coverage += 
								weighted_coverage(multiplex_background_match, multiplex_background_seq);

							// If the multiplex background coverage still looks good, we need to make sure
							// that the *existing* primers do not bind to this proposed amplicon. Up to this
							// point, we have only been checking to ensure that the proposed primers do not
							// bind to the existing amplicons:
							//	  F_a----------R_a <-- Existing
							//  F_b----------------R_b <-- Proposed
							// In the example above, the existing primer F_a and the proposed primer R_b could
							// combine to form a truncated amplicon (that would compete with the intended amplicon
							// F_b + R_b.
							if(s.background_coverage <= opt.max_background_cover){

								// Make a searchable database from the amplicons produced by
								// this assay
								const deque<Sequence> amplicons = 
									trial_assays[t].collect_unique_amplicons(target_keys, target_db, target_seq, 
										opt.target_threshold, opt.target_amplicon_range);

								BitSet local_multiplex_match(amplicons.size(), false);

								for(deque<PCR>::const_iterator i = assay_pool.begin();i != assay_pool.end();++i){

									// Accumulate all of the matches in the local_multiplex_match bitset
									i->find_multiplex_background_match(local_multiplex_match, 
										amplicons, opt, vout);
								}

								// Include the set union of all potential matches between existing primers
								// and the proposed amplicons as background coverage
								s.background_coverage += weighted_coverage(local_multiplex_match, amplicons);
							}
						}
					}

					// This is a candidate for the best assay. Perform detailed screening if there are background
					// sequences present. Note that we first test for background coverage, since the gap-free coverage 
					// computed during the optimization step is a lower bound to the true background coverage. 
					//
					// Important: Since NO GAPS ARE ALLOWED IN TARGET SEQUENCES, we will not be performing a detailed
					// alignment of the assay against the targets (the word-based match in the optimization routine is
	 				// sufficient for targets).
					if(num_active_background > 0){ // There *are* background sequences to test against

						if( (best_score < s) &&
							(s.background_coverage <= opt.max_background_cover) ){

							// Find Smith-Waterman alignment-based amplicon matches
							BitSet background_match(num_background_seq, false);

							// Only match the background sequences. Targets are NOT allowed to have gaps, so our gap free word
							// matching is good enough for counting target matches. The actual identity of target matches will
							// be computed as the end.
							trial_assays[t].find_background_match(background_match, background_keys, background_db,
								background_seq, opt, vout);

							s.background_coverage += weighted_coverage(background_match, background_seq);

							// We already tested for target coverage -- no need to test again
							const bool update_best = (s.background_coverage <= opt.max_background_cover) &&
								(
									(best_score < s) || 
									( (best_score == s) && ( best_assay.total_degeneracy() > trial_assays[t].total_degeneracy() ) )
								);

							#pragma omp critical
							if(update_best){

								best_score = s;
								best_assay.copy_oligos(trial_assays[t]);
								best_background_match = background_match;
							}
						}
					}
					else{ // There are *NO* background sequences to test against
						
						// Note that we still test the value of the background coverage, since multiplex
						// -compatible PCR assays essentially "make their own" background.
						const bool update_best = (s.background_coverage <= opt.max_background_cover) &&
								(
									(best_score < s) || 
									( ( (best_score == s) && ( best_assay.total_degeneracy() > trial_assays[t].total_degeneracy() ) ) )
								);

						#pragma omp critical
						if(update_best){

							best_score = s;
							best_assay.copy_oligos(trial_assays[t]);
						}
					}

					#pragma omp critical
					if(opt.output_filter > Options::VERBOSE){

						vout << "\tCurr accuracy = " << s.accuracy() 
							<< " (" << s.target_coverage << " target, ~" 
							<< s.background_coverage << " background)";

						if(opt.use_multiplex){
							vout << ":" << s.oligo_overlap;
						}

						vout << endl;

						vout << "\tBest accuracy = " << best_score.accuracy() 
							<< " (" << best_score.target_coverage << " target, " 
							<< best_score.background_coverage << " background)";

						if(opt.use_multiplex){
							vout << ": multiplex overlap = " << s.oligo_overlap;
						}

						vout << endl;
					}
				}
			}
						
			// Find the target sequences matched by the best assay. Since gaps between an assay and a target
			// sequence are not allowed, use simple word based matching (as is used in the optimize function).
			BitSet best_target_match;
			deque<Sequence> amplicons;
			deque<AmpliconBounds> bounds; // Amplicon bounds

			// Did we find a "best assay"?
			if(best_score.target_coverage > 0){
			
				best_assay.find_target_match(best_target_match, target_keys, target_db, target_seq, opt);
			
				if(opt.use_multiplex){

					// When designing multiplex PCR-based assays, we need to store the amplicon sequences produced
					// by the best assays. This is to ensure that subsequent primer design do not bind to existing
					// amplicons and potentially inhibit intended PCR.
					//
					// Collect this information now, because each rank has a potentially different
					// set of target words (and only the rank that finds the best assay can actually
					// generate the amplicons for this assay).
					//
					// The bounds deque stores the index and endpoints of each generated amplicon.
					// It will be used to add Base::EOS (i.e. end-of-sequence) symbols to ensure that
					// we do not generate a future assay that flanks a previously designed assay:
					//				F0---------------R0
					//    F1------------------------------R1
					// This flanking scenario would not be avoided simply by adding the first amplicon to
					// the background database, since the new primers (F1 and R1) do not bind to the
					// first amplicon.

					amplicons = 
						best_assay.collect_unique_amplicons(target_keys, target_db, target_seq, 
							opt.target_threshold, opt.target_amplicon_range, &bounds);
				}
			}
			
			// Find the best scoring assay across all MPI ranks
			reduce_best_assay(best_assay, best_score, best_target_match, best_background_match, 
				amplicons, bounds);

			if(best_score.target_coverage <= 0){
				
				// Don't bother writing out assays that failed to detect a single target
				break;
			}

			vout << "\tBest assay: ";
			best_assay.write(vout);
			vout << endl;

			vout << "\tBest accuracy = " << best_score.accuracy()
				<< " (" << best_score.target_coverage << " target, " 
				<< best_score.background_coverage << " background)";

			if(opt.use_multiplex){
				vout << "; multiplex overlap = " << best_score.oligo_overlap;
			}

			vout << endl;

			// Write this assay to disk
			switch(opt.output_format){
				case Options::TEXT_OUTPUT:
					
					fout << "# Assay " << major_assay_id << '.' << minor_assay_id
						<< " has target coverage score = " << best_score.target_coverage 
						<< " (" << (best_score.target_coverage*100.0f)/active_target_norm << "% of active) and background coverage score = " 
						<< best_score.background_coverage << " (" 
						<< ( (num_active_background == 0) ? 0.0f : (best_score.background_coverage*100.0f)/active_background_norm )
						<< "% of active)" << endl;

					fout << "ASSAY." << major_assay_id << '.' << minor_assay_id << '\t';

					break;
				
				case Options::JSON_OUTPUT:
					break;
				
				default:
					throw __FILE__ ":main: Unknown output format";
			};
			
			if(opt.use_multiplex){

				switch(opt.output_format){
					case Options::TEXT_OUTPUT:
					
						best_assay.write(fout, assay_pool);
						fout << endl;
						break;

					case Options::JSON_OUTPUT:
						
						best_assay.write_json(fout, assay_pool);
						break;

					default:
						throw __FILE__ ":main: Unknown output format";
				};

				for(deque<Sequence>::const_iterator i = amplicons.begin();i != amplicons.end();++i){

					i->pack( multiplex_background_db, 
						multiplex_background_seq.size(), /*Used as an index*/
						opt.pack_max_degen, 
						0.0, 1.0, // Don't G+C filter the multiplex background sequences
						opt.min_oligo_length() );

					multiplex_background_seq.push_back(*i);
				}

				// Recompute the multiplex background keys after every new set of amplicons is added
				multiplex_background_keys = keys(multiplex_background_db);

				vout << "\tAdded " << amplicons.size() << " amplicon(s) to the multiplex assay background" << endl;

				// Use the bounds (i.e. the endpoint of each amplicon generated by the current best assay)
				// to add a Base::EOS at the midpoint of each amplicon. This is intended to prevent the selection
				// of a pair of primers that flank the current best assay during a future round of assay design.
				for(deque<AmpliconBounds>::const_iterator i = bounds.begin();i != bounds.end();++i){

					target_seq[i->index].split_sequence(i->begin);

					// Since the amplicon sequences that we put into the multiplex background have their
					// primers trimmed off, we need to also place a split in the center of the amplicon.
					target_seq[i->index].split_sequence( (i->begin + i->end)/2 );

					target_seq[i->index].split_sequence(i->end);
				}
			}
			else{
				switch(opt.output_format){
					case Options::TEXT_OUTPUT:
					
						best_assay.write(fout);
						fout << endl;
						break;

					case Options::JSON_OUTPUT:
						best_assay.write_json(fout);
						break;

					default:
						throw __FILE__ ":main: Unknown output format";
				};
			}

			switch(opt.output_format){
				case Options::TEXT_OUTPUT:

					// Write the sequences matched and mark the matched targets as inactive
					for(unsigned int i = 0;i < num_target_seq;++i){

						if(best_target_match[i]){
							fout << "T-" << target_seq[i].defline() << endl;
						}
					}

					for(unsigned int i = 0;i < num_background_seq;++i){

						if(best_background_match[i]){
							fout << "B-" << background_seq[i].defline() << endl;
						}
					}
					
					break;

				case Options::JSON_OUTPUT:
					
					fout << "\t\t\t\"target matches\":[\n";
					
					{
						bool first_target = true;
						
						for(unsigned int i = 0;i < num_target_seq;++i){

							if(best_target_match[i]){
							
								if(first_target){
									fout << "\t\t\t\t\"" << target_seq[i].defline() << "\"";
								}
								else{
									fout << ",\n\t\t\t\t\"" << target_seq[i].defline() << "\"";
								}
								
								first_target = false;
							}
						}
						
						fout << "\n\t\t\t],\n";
					}
					
					fout << "\t\t\t\"background matches\":[";
					
					{
						bool first_background = true;
						
						for(unsigned int i = 0;i < num_background_seq;++i){

							if(best_background_match[i]){
							
								if(first_background){
									fout << "\n\t\t\t\t\"" << background_seq[i].defline() << "\"";
								}
								else{
									fout << ",\n\t\t\t\t\"" << background_seq[i].defline() << "\"";
								}
								
								first_background = false;
							}
						}
						
						if(first_background){
							fout << "]\n\t\t}";
						}
						else{
							fout << "\n\t\t\t]\n\t\t}";
						}
					}
					
					break;

				default:
					throw __FILE__ ":main: Unknown output format";
			};

			// Mark the matched targets as inactive
			for(unsigned int i = 0;i < num_target_seq;++i){
				if(best_target_match[i]){
					target_seq[i].active(false);
				}
			}
			
			// Save this assay
			assay_pool.push_back(best_assay);
			pool_background.push_back(best_background_match);

			// If we've already reached the desired number of assays, stop searching
			if(assay_iteration >= opt.num_assay){
				break;
			}
		}
		
		if(opt.output_format == Options::JSON_OUTPUT){
			fout << "\n\t],\n";
		}
		
		// Report the number of targets that are *not* detected
		unsigned int undetected_targets = 0;

		for(unsigned int i = 0;i < num_target_seq;++i){

			if( target_seq[i].active() ){
				++undetected_targets;
			}
		}

		if(undetected_targets == 0){
			vout << "Detected all targets" << endl;
		}
		else{
			vout << "Failed to detect a total of " << undetected_targets << " targets" << endl;
		}
		
		// Find all of that background sequences that *are* detected
		BitSet total_background(num_background_seq, false);

		for(deque<BitSet>::const_iterator i = pool_background.begin();i != pool_background.end();++i){
			total_background |= *i;
		}

		const unsigned int num_background_cross_react = total_background.count();

		vout << "Cross reacted with a total of " << num_background_cross_react << " background sequences" << endl;
		
		switch(opt.output_format){
			
			case Options::TEXT_OUTPUT:
			
				fout << "###########################################################################################" << endl;

				if(undetected_targets == 0){
					fout << "# Detected all targets" << endl;
				}
				else{
					fout << "# Failed to detect a total of " << undetected_targets << " targets" << endl;
				}

				if(undetected_targets > 0){

					fout << "# The following targets were *not* detected by any assay" << endl;
					
					for(unsigned int i = 0;i < num_target_seq;++i){

						if( target_seq[i].active() ){
							fout << "-T-" << target_seq[i].defline() << endl;
						}
					}
				}
				
				fout << "###########################################################################################" << endl;
				fout << "# Cross reacted with a total of " << num_background_cross_react << " background sequences" << endl;


				if( num_background_cross_react > 0 ){

					for(unsigned int i = 0;i < num_background_seq;++i){

						if(total_background[i]){
							fout << "+B-" << background_seq[i].defline() << endl;
						}
					}
				}

				break;
			case Options::JSON_OUTPUT:
				
				fout << "\t\"unmatched targets\":[";
				
				if(undetected_targets == 0){
					fout << "],\n";
				}
				else{
					
					bool first_target = true;
					
					for(unsigned int i = 0;i < num_target_seq;++i){

						if( target_seq[i].active() ){
							
							if(first_target){
								fout << "\t\t\"" << target_seq[i].defline() << "\"";
							}
							else{
								fout << ",\n\t\t\"" << target_seq[i].defline() << "\"";
							}
							
							first_target = false;
						}
					}
					
					fout << "\n\t],\n";
				}
				
				fout << "\t\"total number of background matches\":" << num_background_cross_react << ",\n";
				
				if(num_background_cross_react > 0){
				
					fout << "\t\"background matches\":[\n";
					
					bool first_background = true;
					
					for(unsigned int i = 0;i < num_background_seq;++i){

						if(total_background[i]){
							if(first_background){
								fout << "\t\t\"" << background_seq[i].defline() << "\"";
							}
							else{
								fout << ",\n\t\t\"" << background_seq[i].defline() << "\"";
							}
							
							first_background = false;
						}
					}
				}
				else{
					fout << "\t\"background matches\":[]";
				}
				
				fout << "\n}" << endl;
				
				break;
			default:
				throw __FILE__ ":main: Unknown output format";
		};
				
		profile = time(NULL) - profile;
		vout << "Finished in " << time_to_str(profile) << endl;
	}
	catch(const char *error){
		
		cerr << "Caught the error " << error << endl;
		
		MPI_Finalize();
		
		return EXIT_FAILURE;
	}
	catch(const string error){
		
		cerr << "Caught the error " << error << endl;
		
		MPI_Finalize();
		
		return EXIT_FAILURE;
	}
	catch(...){
	
		cerr << "Caught an unhandled error" << endl;
		
		MPI_Finalize();
		
		return EXIT_FAILURE;
	}
	
	MPI_Finalize();
	
	return EXIT_SUCCESS;
}

string time_to_str(const time_t &m_elapsed)
{
	const float elapsed_sec = m_elapsed % 60;
	
	float elapsed_time = (m_elapsed - elapsed_sec)/60.0f; // In min
	
	const float elapsed_min = fmod(elapsed_time, 60.0f);
	elapsed_time = (elapsed_time - elapsed_min)/60.0f; // In hour
	
	const float elapsed_hour = fmod(elapsed_time, 24.0f);
	elapsed_time = (elapsed_time - elapsed_hour)/24.0f; // In day
	
	stringstream sout;
	
	sout << "Run time is " 
		<< elapsed_time 
		<< " days, "
		<< elapsed_hour
		<< ( (elapsed_hour == 1) ? " hour, " : " hours, ")
		<< elapsed_min
		<< " min and "
		<< elapsed_sec
		<< " sec";
	
	return sout.str();
}

void sequence_summary(const string &m_prefix, const deque<Sequence> &m_seq, ostream &m_out,
	const Options::OutputFormat &m_format)
{
	float ave = 0.0;
	float stdev = 0.0;
	
	switch(m_format){
		case Options::TEXT_OUTPUT:
			m_out << m_prefix <<  " Number of sequences = " << m_seq.size() << endl;
			break;
		case Options::JSON_OUTPUT:
			m_out << "\t\"" << m_prefix << "\":{\n"
				<< "\t\t\"number of sequences\":" << m_seq.size();
			break;
		default:
			throw __FILE__ ":sequence_summary: Unknown output format";
	};
	
	if( m_seq.empty() ){
		
		if(m_format == Options::JSON_OUTPUT){
			m_out << "\n\t},\n";
		}
		
		return;
	}
	
	unsigned int min_length = m_seq[0].length();
	unsigned int max_length = m_seq[0].length();
	
	for(deque<Sequence>::const_iterator i = m_seq.begin();i != m_seq.end();++i){
	
		const unsigned int l = i->length();
		
		ave += l;
		
		min_length = min(min_length, l);
		max_length = max(max_length, l);
	}
	
	ave /= m_seq.size();
	
	for(deque<Sequence>::const_iterator i = m_seq.begin();i != m_seq.end();++i){
		
		const float tmp = i->length() - ave;
		
		stdev += tmp*tmp;
	}

	if( m_seq.size() > 1){	
		stdev = sqrt( stdev/(m_seq.size() - 1) );
	}
	else{
		stdev = 0.0;
	}
	
	switch(m_format){
		case Options::TEXT_OUTPUT:
		
			m_out << m_prefix << " Min sequence length = " << min_length << endl;
			m_out << m_prefix << " Max sequence length = " << max_length << endl;
			m_out << m_prefix << " Average sequence length = " << ave << endl;
			m_out << m_prefix << " Stdev sequence length = " << stdev << endl;
			break;
		case Options::JSON_OUTPUT:
		
			m_out << ",\n\t\t\"min sequence length\":" << min_length << ",\n"
				<< "\t\t\"max sequence length\":" << max_length << ",\n"
				<< "\t\t\"average sequence length\":" << ave << ",\n"
				<< "\t\t\"stdev sequence length\":" << stdev << "\n\t},\n";
			break;
		default:
			throw __FILE__ ":sequence_summary: Unknown output format";
	};
}

float weighted_coverage(const BitSet &m_match, const deque<Sequence> &m_seq)
{
	double ret = 0.0; // Accumulate in double
	
	const size_t len = m_seq.size();
	
	assert( len == m_match.size() );
	
	for(size_t i = 0;i < len;++i){
		
		if(m_match[i]){
			ret += m_seq[i].weight();
		}
	}
	
	return ret;
}

// Return the rank that owns the best assay.
void reduce_best_assay(PCR &m_assay, Score &m_score, BitSet &m_target_match, BitSet &m_background_match,
	deque<Sequence> &m_amplicons, deque<AmpliconBounds> &m_bounds)
{
	if(mpi_numtasks <= 1){
		return;
	}
	
	const bool is_root = (mpi_rank == 0);
				
	// Have every worker send their assay/scores/background match data to the id 0 task
	if(is_root){
		
		MPI_Status status;
		
		for(int task = 1;task < mpi_numtasks;++task){
			
			unsigned int len = 0;
			
			if(MPI_Recv(&len, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, BEST_ASSAY, MPI_COMM_WORLD, &status) != MPI_SUCCESS){
				throw __FILE__ ":reduce_best_assay: Error receiving msg";
			}
			
			unsigned char *buffer = new unsigned char[len];
			
			if(buffer == NULL){
				throw __FILE__ ":reduce_best_assay: Unable to allocate receive buffer";
			}
			
			if(MPI_Recv(buffer, len, MPI_BYTE, status.MPI_SOURCE, BEST_ASSAY_DATA, MPI_COMM_WORLD, &status) != MPI_SUCCESS){
				throw __FILE__ ":reduce_best_assay: Error receiving data";
			}
			
			unsigned char* ptr = buffer;
			
			// Unpack the score
			Score trial_score;
			
			ptr = mpi_unpack(ptr, trial_score);
			
			if(trial_score < m_score){
				
				delete [] buffer;
				continue;
			}
			
			// Unpack the assay degeneracy
			double trial_degeneracy;
			
			memcpy( &trial_degeneracy, ptr, sizeof(trial_degeneracy) );
			ptr += sizeof(trial_degeneracy);
			
			// In the case of a tied score, save the lowest degeneracy assay as the best
			if( (trial_score == m_score) && (m_assay.total_degeneracy() <= trial_degeneracy) ){
				
				delete [] buffer;
				continue;
			}
			
			m_score = trial_score;
			
			// Unpack the assay and the target & background bitsets
			ptr = m_assay.mpi_unpack(ptr);
			
			ptr = mpi_unpack(ptr, m_target_match);
			
			ptr = mpi_unpack(ptr, m_background_match);
			
			ptr = mpi_unpack(ptr, m_amplicons);

			ptr = mpi_unpack(ptr, m_bounds);
			
			delete [] buffer;
		}
		
		////////////////////////////////////////////////////////////////
		// The root task now sends out the best assay data to all tasks
		////////////////////////////////////////////////////////////////
		
		unsigned int len = mpi_size(m_score) + m_assay.mpi_size() + 
			mpi_size(m_target_match) + mpi_size(m_background_match) +
			mpi_size(m_amplicons) + mpi_size(m_bounds);
		
		MPI_Bcast( (void*)&(len), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
		
		unsigned char* buffer = new unsigned char [len];
		
		if(buffer == NULL){
			throw __FILE__ ":reduce_best_assay: Unable to allocate send buffer (0)";
		}
		
		unsigned char *ptr = buffer;
		
		ptr = mpi_pack(ptr, m_score);
		
		ptr = m_assay.mpi_pack(ptr);
		
		ptr = mpi_pack(ptr, m_target_match);
		
		ptr = mpi_pack(ptr, m_background_match);
		
		ptr = mpi_pack(ptr, m_amplicons);

		ptr = mpi_pack(ptr, m_bounds);
		
		MPI_Bcast(buffer, len, MPI_BYTE, 0, MPI_COMM_WORLD);
		
		delete [] buffer;
	}
	else{
		
		// Pack a buffer with all of the data to send to id 0
		unsigned int len = mpi_size(m_score) + sizeof(double) + m_assay.mpi_size() + 
			mpi_size(m_target_match) + mpi_size(m_background_match) + 
			mpi_size(m_amplicons) + mpi_size(m_bounds);
		
		unsigned char* buffer = new unsigned char [len];
		
		if(buffer == NULL){
			throw __FILE__ ":reduce_best_assay: Unable to allocate send buffer";
		}
		
		unsigned char *ptr = buffer;
		
		ptr = mpi_pack(ptr, m_score);
		
		double local_degeneracy = m_assay.total_degeneracy();
		
		memcpy( ptr, &local_degeneracy, sizeof(local_degeneracy) );
		ptr += sizeof(local_degeneracy);
		
		ptr = m_assay.mpi_pack(ptr);
		
		ptr = mpi_pack(ptr, m_target_match);
		
		ptr = mpi_pack(ptr, m_background_match);
		
		ptr = mpi_pack(ptr, m_amplicons);

		ptr = mpi_pack(ptr, m_bounds);
		
		if(MPI_Send( (void*)&len, 1, MPI_UNSIGNED, 0, BEST_ASSAY, MPI_COMM_WORLD ) != MPI_SUCCESS){
			throw __FILE__ ":reduce_best_assay: Error sending msg";
		}
		
		if(MPI_Send(buffer, len, MPI_BYTE, 0, BEST_ASSAY_DATA, MPI_COMM_WORLD ) != MPI_SUCCESS){
			throw __FILE__ ":reduce_best_assay: Error sending data";
		}
		
		delete [] buffer;
		
		///////////////////////////////////////////////////////////////////////
		// All non-root tasks now receive the best assay data from the root
		///////////////////////////////////////////////////////////////////////
		
		MPI_Bcast( (void*)&(len), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD );
		
		buffer = new unsigned char [len];
		
		if(buffer == NULL){
			throw __FILE__ ":reduce_best_assay: Unable to allocate send buffer (non-0)";
		}
		
		MPI_Bcast(buffer, len, MPI_BYTE, 0, MPI_COMM_WORLD);
		
		ptr = buffer;
		
		ptr = mpi_unpack(ptr, m_score);
		
		ptr = m_assay.mpi_unpack(ptr);
		
		ptr = mpi_unpack(ptr, m_target_match);
		
		ptr = mpi_unpack(ptr, m_background_match);
		
		ptr = mpi_unpack(ptr, m_amplicons);
		
		ptr = mpi_unpack(ptr, m_bounds);

		delete [] buffer;
	}
}

string truncate_prefix(const string &m_str, const size_t &m_max_len)
{
	if(m_str.size() <= m_max_len){
		return m_str;
	}
	
	// Is there room for "..."?
	if(m_max_len <= 3){
		return m_str.substr(m_str.size() - m_max_len, m_max_len);
	}
	
	return string("...") + m_str.substr( (m_str.size() + 3) - m_max_len, m_max_len - 3);
}
