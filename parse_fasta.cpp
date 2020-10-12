#include "pcramp.h"
#include "update.h"
#include <fstream>

#include <zlib.h>

using namespace std;

void parse_fasta(const string &m_filename, deque<Sequence> &m_dbase,
	const size_t &m_min_length, const size_t &m_max_length,
	const deque<string> &m_ignore)
{
	
	// Use zlib to read both compressed and uncompressed fasta files.
	gzFile fin = gzopen(m_filename.c_str(), "r");
	
	if(fin == NULL){
		
		cerr << "Error opening: " << m_filename << endl;
		throw __FILE__ ":parse_fasta: Unable to open fasta file";
	}
	
	const int buffer_len = 2048;
	char buffer[buffer_len];
	
	string defline;
	deque<char> seq;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		if(strchr(buffer, '>') != NULL){
			
			if( !seq.empty() ){

				if( (seq.size() >= m_min_length) && (seq.size() <= m_max_length) &&
				    !ignore_record(defline, m_ignore) ){
				
					m_dbase.push_back( Sequence() );

					Sequence &ref = m_dbase.back();

					ref = seq; // Only copies the sequences data
					ref.defline(defline); // Sets the weight from the defline (if possible)
					ref.active(true);
				}
				
				defline.clear();
				seq.clear();
			}
			
			// Remove any end of line symbols
			for(char* p = buffer;*p != '\0';++p){
				if( (*p == '\n') || (*p == '\r') ){
					*p = '\0';
				}
			}

			defline = buffer;
		}
		else{
			for(char* p = buffer;*p != '\0';++p){
				
				if( !isspace(*p) ){
				
					// Preserve the case of all sequences
					// (since upper and lower case bases may receive
					// different weights)
					seq.push_back(*p);
				}
			}
		}
	}
	
	if( (seq.size() >= m_min_length) && (seq.size() <= m_max_length) &&
	    !ignore_record(defline, m_ignore) ){

		m_dbase.push_back( Sequence() );

		Sequence &ref = m_dbase.back();

		ref = seq; // Only copies the sequences data
		ref.defline(defline); // Sets the weight from the defline (if possible)
		ref.active(true);
	}
	
	gzclose(fin);
}

void append_fasta_group(const string &m_filename, Sequence &m_seq,
	const size_t &m_min_length, const size_t &m_max_length,
	const size_t &m_num_pad, const deque<string> &m_ignore)
{
	
	// Use zlib to read both compressed and uncompressed fasta files.
	gzFile fin = gzopen(m_filename.c_str(), "r");
	
	if(fin == NULL){
		
		cerr << "Error opening: " << m_filename << endl;
		throw __FILE__ ":append_fasta_group: Unable to open fasta file";
	}
	
	const int buffer_len = 2048;
	char buffer[buffer_len];
	
	string defline;
	deque<char> seq;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		if(strchr(buffer, '>') != NULL){
			
			if( !seq.empty() ){

				if( (seq.size() >= m_min_length) && (seq.size() <= m_max_length) &&
				    !ignore_record(defline, m_ignore) ){
				
					if( !m_seq.empty() ){
						
						// Append m_num_pad Base::EOS characters between 
						// each fasta record
						m_seq.pad(m_num_pad);
					}
					
					m_seq += seq; // Append to the current sequence
				}
				
				defline.clear();
				seq.clear();
			}
			
			// Remove any end of line symbols
			for(char* p = buffer;*p != '\0';++p){
				if( (*p == '\n') || (*p == '\r') ){
					*p = '\0';
				}
			}

			defline = buffer;
		}
		else{
			for(char* p = buffer;*p != '\0';++p){
				
				if( !isspace(*p) ){
				
					// Preserve the case of all sequences
					// (since upper and lower case bases may receive
					// different weights)
					seq.push_back(*p);
				}
			}
		}
	}
	
	if( (seq.size() >= m_min_length) && (seq.size() <= m_max_length) &&
	    !ignore_record(defline, m_ignore) ){
	    
	    	if( !m_seq.empty() ){
						
			// Append m_num_pad Base::EOS charaters between 
			// each fasta record
			m_seq.pad(m_num_pad);
		}
		
		m_seq += seq; // Append to the current sequence
	}
	
	gzclose(fin);
}

bool ignore_record(const string &m_defline, const deque<string> &m_ignore)
{
	if( m_ignore.empty() ){
		return false;
	}
	
	// The strings in the m_ignore structure are already lower case.
	const string def = tolower(m_defline);
	
	for(deque<string>::const_iterator i = m_ignore.begin();i != m_ignore.end();++i){
		
		if( def.find(*i) != string::npos){
			return true;
		}
	}
	
	return false;
}
