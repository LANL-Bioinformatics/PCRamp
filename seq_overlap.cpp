#include "seq_overlap.h"
#include <math.h>
#include <iostream>
#include <sstream>
#include <algorithm>

#ifdef DETECT_CPU
#include <cpuid.h>
#endif // DETECT_CPU

using namespace std;
using namespace SO;

SeqOverlap::SeqOverlap(const AlignmentMode &m_mode, const bool &m_is_na) :
	is_na(m_is_na)
{	

	#ifdef DETECT_CPU
	// Use the cpuid function (macro) to check the capabilities
	// of the cpu we're running on
	unsigned int level = 1;
	
	unsigned int eax;
	unsigned int ebx;
	unsigned int ecx;
	unsigned int edx;
	
	const int ret = __get_cpuid(level, &eax, &ebx, &ecx, &edx);
	
	cerr << "cpuid ret = " << ret << endl;
	cerr << "eax = " << eax << endl;
	cerr << "ebx = " << ebx << endl;
	cerr << "ecx = " << ecx << endl;
	cerr << "edx = " << edx << endl;
	
	cerr << "bit_MMX = " << (bit_MMX & edx ? "true" : "false") << endl;
	cerr << "bit_SSE = " << (bit_SSE & edx ? "true" : "false") << endl;
	cerr << "bit_SSE2 = " << (bit_MMX & edx ? "true" : "false") << endl;
	cerr << "bit_SSE3 = " << (bit_SSE3 & ecx  ? "true" : "false") << endl;
	cerr << "bit_SSE4a = " << (bit_SSE4a & ecx  ? "true" : "false") << endl;
	cerr << "bit_SSE4_2 = " << (bit_SSE4_2 & ecx ? "true" : "false") << endl;
	cerr << "bit_SSE5 = " << (bit_SSE5 & ecx ? "true" : "false") << endl;
	cerr << "bit_AVX = " << (bit_AVX & ecx ? "true" : "false") << endl;
	#endif // DETECT_CPU
	
	mode = m_mode;
	
	mask_N_na = false;
	
	last_row = curr_row = NULL;
	dp_num_elem = 0;
	
	max_query_len = 0;
	max_target_len = 0;
	
	query_len.sse = _MM_SET1(0);
	target_len.sse = _MM_SET1(0);
	
	if(is_na){
	
		// These are blastn defaults, note that other nucleotide alignment algorithms (like FASTA, megablast, etc.) 
		// use different paramters
		match.sse = _MM_SET1(2);
		mask.sse = _MM_SET1(0);
		mismatch.sse = _MM_SET1(-3);
		gap_existance.sse = _MM_SET1(-5);
		gap_extension.sse = _MM_SET1(-2);
	}
	else{
		// Gaps costs from BLASTP default parameters
		gap_existance.sse = _MM_SET1(-11);
		gap_extension.sse = _MM_SET1(-1);
		
		// There are number of possible protein score matricies: 
		// PAM30, PAM70,  BLOSUM45, BLOSUM62, BLOSUM80
		// By default, BLASTP used BLOSUM62
		init_BLOSUM62(aa_score_matrix);
		
		//throw __FILE__ ": Amino acid alignments are not currently enabled!";
	}
	
	// Allow dynamic sizing for both the target and query
	query = (sse_elem*)_mm_malloc(DEFAULT_SEQUENCE_LENGTH*sizeof(sse_elem), 16);
	
	if(query == NULL){
		throw __FILE__ ":SeqOverlap::SeqOverlap: Unable to allocate query sequence storage";
	}
	
	target = (sse_elem*)_mm_malloc(DEFAULT_SEQUENCE_LENGTH*sizeof(sse_elem), 16);
	
	if(target == NULL){
		throw __FILE__ ":SeqOverlap::SeqOverlap: Unable to allocate target sequence storage";
	}
	
	max_allowed_query_len = DEFAULT_SEQUENCE_LENGTH;
	max_allowed_target_len = DEFAULT_SEQUENCE_LENGTH;
};

// Compute the alignment between two sequences; the query 
// and the target. Both the query and target sequences are assumed to be
// in 5'-3' orientation.
void SeqOverlap::align_overlap()
{
	// Add one to allow for an initial row and column of zeros
	const size_t num_elem = max_target_len + 1;

	const sse_elem zero(0);
	const sse_elem all_N(NA::N);
	const sse_elem m_minus_mm( _MM_SUB(match.sse, mismatch.sse) );
	const sse_elem query_length_minus_1( _MM_SUB( query_len.sse, _MM_SET1(1) ) );
	const sse_elem target_length_minus_1( _MM_SUB( target_len.sse, _MM_SET1(1) ) );
	
	// Resize the dynamic programming matrix if needed
	if( !last_row || !curr_row || (num_elem > dp_num_elem) ){

		if(last_row){
			_mm_free(last_row);
		}
		
		if(curr_row){
			_mm_free(curr_row);
		}
		
		dp_num_elem = num_elem;
		
		last_row = (SO_Elem*)_mm_malloc(dp_num_elem*sizeof(SO_Elem), 16);
		
		if(!last_row){
			throw __FILE__ ":SeqOverlap::align_overlap: Unable to allocate memory for last_row";
		}

		curr_row = (SO_Elem*)_mm_malloc(dp_num_elem*sizeof(SO_Elem), 16);
		
		if(!curr_row){
			throw __FILE__ ":SeqOverlap::align_overlap: Unable to allocate memory for curr_row";
		}
	}
	
	// Initialize the dynamic programming matrix to have zeros along the first row 
	for(SO_Score j = 0;j <= max_target_len;j++){

		SO_Elem &elem_ref = last_row[j];

		elem_ref.M.sse = zero.sse;
		elem_ref.I_query.sse = elem_ref.I_target.sse = gap_existance.sse;
		
		elem_ref.M_start_i.sse = zero.sse;
		elem_ref.M_start_j.sse = _MM_SET1(j);		
	}
	
	// Prevent valgrind "Conditional jump or move depends on uninitialised value(s)" error.
	// I belive that this error is spurious, but I'm not 100% sure ...
	memset( curr_row, 0x0, dp_num_elem*sizeof(SO_Elem) );
	
	// Reset the maximum element
	max_elem.M.sse = _MM_SET1(MINIMUM_ALLOWED_SCORE);
	
	for(SO_Score i = 0;i < max_query_len;i++){
		
		// Initialize the dynamic programming matrix to have zeros along the first column
		SO_Elem &elem_ref = curr_row[0];

		elem_ref.M.sse = zero.sse;
		elem_ref.I_query.sse = elem_ref.I_target.sse = gap_existance.sse;
		
		elem_ref.M_start_i.sse = _MM_SET1(i + 1);
		elem_ref.M_start_j.sse = zero.sse;
		
		// The dp matrix has max_query_len rows and max_target_len columns
		// A B
		// C X <-- dp[i][j]
		SO_Elem *A_ptr = last_row;
		SO_Elem *B_ptr = A_ptr + 1;
		SO_Elem *C_ptr = curr_row;
		SO_Elem *X_ptr = C_ptr + 1;
		
		const sse_elem q = query[i];
		
		const sse_elem all_i( _MM_SET1(i) );
		const sse_elem valid_query( _MM_CMPLT(all_i.sse, query_len.sse) );
		const sse_elem query_is_N( _MM_CMPEQ(q.sse, all_N.sse) );
		 
		for(SO_Score j = 0;j < max_target_len;j++, A_ptr++, B_ptr++, C_ptr++, X_ptr++){
			
			sse_elem tmp_a;
			sse_elem tmp_b;
			sse_elem tmp_c;
			
			const sse_elem t = target[j];
			
			if(is_na){
			
				tmp_a.sse = _MM_ADD(
						_mm_and_si128(
							_MM_CMPGT(_mm_and_si128(q.sse, target[j].sse), zero.sse),
							m_minus_mm.sse),
						mismatch.sse);

				if(mask_N_na){

					tmp_b.sse = _mm_or_si128(query_is_N.sse, _MM_CMPEQ(target[j].sse, all_N.sse) );

					tmp_a.sse = _mm_or_si128( 
							_mm_and_si128(tmp_b.sse, mask.sse), 
							_mm_andnot_si128(tmp_b.sse, tmp_a.sse) );
				}

				X_ptr->M.sse = _MM_ADD(
					_MM_MAX( _MM_MAX(A_ptr->M.sse, A_ptr->I_query.sse), A_ptr->I_target.sse ),
					 tmp_a.sse);
			}
			else{
				
				// See the align_smith_waterman() for comments on this serial gather
				for(unsigned char local = 0;local < SO_LEN;local++){
					X_ptr->M.v[local] = aa_score_matrix[ q.v[local]*AA::NUM_AA + t.v[local] ];
				}
				
				X_ptr->M.sse = _MM_ADD(
					_MM_MAX( _MM_MAX(A_ptr->M.sse, A_ptr->I_query.sse), A_ptr->I_target.sse ),
					X_ptr->M.sse);
			}
			
			//////////////////////////////////////////////////////////////////////////////
			tmp_a.sse = _mm_or_si128(
					_MM_CMPLT(A_ptr->M.sse, A_ptr->I_query.sse), 
					_MM_CMPLT(A_ptr->M.sse, A_ptr->I_target.sse) );
					
			X_ptr->M_start_i.sse = _mm_andnot_si128(tmp_a.sse, A_ptr->M_start_i.sse);
			X_ptr->M_start_j.sse = _mm_andnot_si128(tmp_a.sse, A_ptr->M_start_j.sse);
			
			//////////////////////////////////////////////////////////////////////////////
			tmp_a.sse = _mm_andnot_si128(
					_MM_CMPLT(A_ptr->I_query.sse, A_ptr->I_target.sse),
					_MM_CMPGT(A_ptr->I_query.sse, A_ptr->M.sse) );
			
			X_ptr->M_start_i.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_i.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_query_start_i.sse) );
			
			X_ptr->M_start_j.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_j.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_query_start_j.sse) );
			
			//////////////////////////////////////////////////////////////////////////////
			tmp_a.sse = _mm_and_si128(
					_MM_CMPGT(A_ptr->I_target.sse, A_ptr->M.sse),
					_MM_CMPGT(A_ptr->I_target.sse, A_ptr->I_query.sse) );
						
			X_ptr->M_start_i.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_i.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_target_start_i.sse) );
			
			X_ptr->M_start_j.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_j.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_target_start_j.sse) );
			
			//////////////////////////////////////////////////////////////////////////////////////
			// Unlike the SmithWaterman alignment, *don't* clamp the C cell values to zero
			tmp_b.sse = _MM_ADD(C_ptr->M.sse, gap_existance.sse);
			tmp_c.sse = _MM_ADD(C_ptr->I_query.sse, gap_extension.sse);
			
			X_ptr->I_query.sse = _MM_MAX(tmp_b.sse, tmp_c.sse );

			tmp_a.sse = _MM_CMPLT(tmp_b.sse, tmp_c.sse);
			
			X_ptr->I_query_start_i.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, C_ptr->M_start_i.sse),
							_mm_and_si128(
								tmp_a.sse, C_ptr->I_query_start_i.sse) );
			
			X_ptr->I_query_start_j.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, C_ptr->M_start_j.sse),
							_mm_and_si128(
								tmp_a.sse, C_ptr->I_query_start_j.sse) );
			
			//////////////////////////////////////////////////////////////////////////////////////
			// Unlike the SmithWaterman alignment, *don't* clamp the B cell values to zero
			tmp_b.sse = _MM_ADD(B_ptr->M.sse, gap_existance.sse);
			tmp_c.sse = _MM_ADD(B_ptr->I_target.sse, gap_extension.sse);
								
			X_ptr->I_target.sse = _MM_MAX(tmp_b.sse, tmp_c.sse);
			
			tmp_a.sse = _MM_CMPLT(tmp_b.sse, tmp_c.sse);
			
			X_ptr->I_target_start_i.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, B_ptr->M_start_i.sse),
							_mm_and_si128(
								tmp_a.sse, B_ptr->I_target_start_i.sse) );
			
			X_ptr->I_target_start_j.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, B_ptr->M_start_j.sse),
							_mm_and_si128(
								tmp_a.sse, B_ptr->I_target_start_j.sse) );
								
			//////////////////////////////////////////////////////////////////////////////
			
			const sse_elem all_j( _MM_SET1(j) );
			const sse_elem valid_query_and_target( 
				_mm_and_si128(valid_query.sse, _MM_CMPLT(all_j.sse, target_len.sse) ) );
				
			//  ((i == query size - 1) || (j == target_size - 1)) && (X_ptr->M >= max_elem.M)
			tmp_a.sse = _mm_and_si128(valid_query_and_target.sse,
					_mm_andnot_si128(
						_MM_CMPLT(X_ptr->M.sse, max_elem.M.sse),
						_mm_or_si128(
							_MM_CMPEQ(all_i.sse, query_length_minus_1.sse),
							_MM_CMPEQ(all_j.sse, target_length_minus_1.sse) ) ) );
			
			max_elem.M.sse = _mm_or_si128(
						_mm_andnot_si128(tmp_a.sse, max_elem.M.sse),
						_mm_and_si128(tmp_a.sse, X_ptr->M.sse) );
			
			max_elem.M_start_i.sse = _mm_or_si128(
						_mm_andnot_si128(tmp_a.sse, max_elem.M_start_i.sse),
						_mm_and_si128(tmp_a.sse, X_ptr->M_start_i.sse) );
			
			max_elem.M_start_j.sse = _mm_or_si128(
						_mm_andnot_si128(tmp_a.sse, max_elem.M_start_j.sse),
						_mm_and_si128(tmp_a.sse, X_ptr->M_start_j.sse) );
							
			stop_i.sse = _mm_or_si128(
						_mm_andnot_si128(tmp_a.sse, stop_i.sse),
						_mm_and_si128(tmp_a.sse, all_i.sse) );
			
			stop_j.sse = _mm_or_si128(
						_mm_andnot_si128(tmp_a.sse, stop_j.sse),
						_mm_and_si128(tmp_a.sse, all_j.sse) );
		}
		
		// Swap the last_row and curr_row
		swap(last_row, curr_row);
	}
}

// Compute the alignment between two sequences; the query 
// and the target. Both the query and target sequences are assumed to be
// in 5'-3' orientation.
void SeqOverlap::align_smith_waterman()
{
	// Add one to allow for an initial row and column of zeros
	const size_t num_elem = max_target_len + 1;

	const sse_elem zero(0);
	const sse_elem all_N(NA::N);
	const sse_elem m_minus_mm( _MM_SUB(match.sse, mismatch.sse) );
	
	// Resize the dynamic programming matrix if needed
	if( !last_row || !curr_row || (num_elem > dp_num_elem) ){

		if(last_row){
			_mm_free(last_row);
		}
		
		if(curr_row){
			_mm_free(curr_row);
		}
		
		dp_num_elem = num_elem;
		
		last_row = (SO_Elem*)_mm_malloc(dp_num_elem*sizeof(SO_Elem), 16);
		
		if(!last_row){
			throw __FILE__ ":SeqOverlap::align_overlap: Unable to allocate memory for last_row";
		}
		
		curr_row = (SO_Elem*)_mm_malloc(dp_num_elem*sizeof(SO_Elem), 16);
		
		if(!curr_row){
			throw __FILE__ ":SeqOverlap::align_overlap: Unable to allocate memory for curr_row";
		}
	}

	// Initialize the dynamic programming matrix to have zeros along the first row 
	for(SO_Score j = 0;j <= max_target_len;j++){

		SO_Elem &elem_ref = last_row[j];

		elem_ref.M.sse = zero.sse;
		elem_ref.I_query.sse = elem_ref.I_target.sse = gap_existance.sse;
		
		elem_ref.M_start_i.sse = zero.sse;
		elem_ref.M_start_j.sse = _MM_SET1(j);		
	}
	
	// Prevent valgrind "Conditional jump or move depends on uninitialised value(s)" error.
	// I belive that this error is spurious, but I'm not 100% sure ...
	memset( curr_row, 0x0, dp_num_elem*sizeof(SO_Elem) );
	
	// Reset the maximum element
	max_elem.M.sse = _MM_SET1(0);
	
	for(SO_Score i = 0;i < max_query_len;i++){
		
		// Initialize the dynamic programming matrix to have zeros along the first column
		SO_Elem &elem_ref = curr_row[0];

		elem_ref.M.sse = zero.sse;
		elem_ref.I_query.sse = elem_ref.I_target.sse = gap_existance.sse;
		
		//elem_ref.M_start = make_pair(i + 1, 0);
		elem_ref.M_start_i.sse = _MM_SET1(i + 1);
		elem_ref.M_start_j.sse = zero.sse;
		
		// The dp matrix has max_query_len rows and max_target_len columns
		// A B
		// C X <-- dp[i][j]
		SO_Elem *A_ptr = last_row;
		SO_Elem *B_ptr = A_ptr + 1;
		SO_Elem *C_ptr = curr_row;
		SO_Elem *X_ptr = C_ptr + 1;
		
		const sse_elem q = query[i];
		
		const sse_elem all_i(i);
		const sse_elem valid_query( _MM_CMPLT(all_i.sse, query_len.sse) );
		const sse_elem query_is_N( _MM_CMPEQ(q.sse, all_N.sse) );
		
		for(SO_Score j = 0;j < max_target_len;j++, A_ptr++, B_ptr++, C_ptr++, X_ptr++){
			
			sse_elem tmp_a;
			sse_elem tmp_b;
			sse_elem tmp_c;
			
			const sse_elem all_j(j);
			const sse_elem t = target[j];
			
			tmp_c.sse = _MM_MAX( _MM_MAX(A_ptr->M.sse, A_ptr->I_query.sse), A_ptr->I_target.sse );
			
			if(is_na){
			
				tmp_a.sse = _MM_ADD(
					_mm_and_si128(
						_MM_CMPGT(_mm_and_si128(q.sse, t.sse), zero.sse),
						m_minus_mm.sse),
					mismatch.sse);
					
				if(mask_N_na){
				
					tmp_b.sse = _mm_or_si128(query_is_N.sse, _MM_CMPEQ(t.sse, all_N.sse) );

					tmp_a.sse = _mm_or_si128( 
							_mm_and_si128(tmp_b.sse, mask.sse), 
							_mm_andnot_si128(tmp_b.sse, tmp_a.sse) );
				}
				
				X_ptr->M.sse = _MM_ADD(
					// Unlike the Overlap alignment, clamp the max value at zero
					_MM_MAX(tmp_c.sse, zero.sse),
					tmp_a.sse);
			}
			else{
				// Interesting! The naive serial version of gathering scores from the amino acid scoring
				// matrix appears to be faster than the slightly parallel version
				// Totally serial version of loading from the amino acid scoring matrix (parallel speed up ~ 5.4).
				// Note that some future version of AVX will likely support SIMD gather
				for(unsigned char local = 0;local < SO_LEN;local++){
					X_ptr->M.v[local] = aa_score_matrix[ q.v[local]*AA::NUM_AA + t.v[local] ];
				}
				
				// Partially parallel version (parallel speed up ~ 4.2)
				//X_ptr->M.sse = _MM_ADD(_MM_MULLO(q.sse, all_NUM_AA.sse), t.sse);
				
				//for(unsigned char local = 0;local < SO_LEN;local++){
				//	X_ptr->M.v[local] = aa_score_matrix[ X_ptr->M.v[local] ];
				//}
				
				X_ptr->M.sse = _MM_ADD(
					// Unlike the Overlap alignment, clamp the max value at zero
					_MM_MAX(tmp_c.sse, zero.sse),
					X_ptr->M.sse);
			}
			
			//////////////////////////////////////////////////////////////////////////////
			tmp_a.sse = _mm_or_si128(
					_MM_CMPLT(A_ptr->M.sse, A_ptr->I_query.sse), 
					_MM_CMPLT(A_ptr->M.sse, A_ptr->I_target.sse) );
								
			X_ptr->M_start_i.sse = _mm_andnot_si128(tmp_a.sse, A_ptr->M_start_i.sse);
			X_ptr->M_start_j.sse = _mm_andnot_si128(tmp_a.sse, A_ptr->M_start_j.sse);
			
			//////////////////////////////////////////////////////////////////////////////
			tmp_a.sse = _mm_andnot_si128(
					_MM_CMPLT(A_ptr->I_query.sse, A_ptr->I_target.sse),
					_MM_CMPGT(A_ptr->I_query.sse, A_ptr->M.sse) );
			
			X_ptr->M_start_i.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_i.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_query_start_i.sse) );
			
			X_ptr->M_start_j.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_j.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_query_start_j.sse) );
			
			//////////////////////////////////////////////////////////////////////////////
			tmp_a.sse = _mm_and_si128(
					_MM_CMPGT(A_ptr->I_target.sse, A_ptr->M.sse),
					_MM_CMPGT(A_ptr->I_target.sse, A_ptr->I_query.sse) );
						
			X_ptr->M_start_i.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_i.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_target_start_i.sse) );
			
			X_ptr->M_start_j.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_j.sse),
							_mm_and_si128(tmp_a.sse, A_ptr->I_target_start_j.sse) );
			
			//////////////////////////////////////////////////////////////////////////////
			// Check for values that were clamped to zero
			tmp_a.sse = _MM_CMPGT(zero.sse, tmp_c.sse);
			
			X_ptr->M_start_i.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_i.sse),
							_mm_and_si128(tmp_a.sse, all_i.sse) );
			
			X_ptr->M_start_j.sse = _mm_or_si128( 
							_mm_andnot_si128(
								tmp_a.sse, X_ptr->M_start_j.sse),
							_mm_and_si128(tmp_a.sse, all_j.sse) );
							
			//////////////////////////////////////////////////////////////////////////////
			// Unlike the Overlap alignment, clamp the C cell values to zero
			tmp_b.sse = _MM_ADD(_MM_MAX(C_ptr->M.sse, zero.sse), gap_existance.sse);
			tmp_c.sse = _MM_ADD(_MM_MAX(C_ptr->I_query.sse, zero.sse), gap_extension.sse);
			
			X_ptr->I_query.sse = _MM_MAX(tmp_b.sse, tmp_c.sse );

			tmp_a.sse = _MM_CMPLT(tmp_b.sse, tmp_c.sse);
			
			X_ptr->I_query_start_i.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, C_ptr->M_start_i.sse),
							_mm_and_si128(
								tmp_a.sse, C_ptr->I_query_start_i.sse) );
			
			X_ptr->I_query_start_j.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, C_ptr->M_start_j.sse),
							_mm_and_si128(
								tmp_a.sse, C_ptr->I_query_start_j.sse) );
			
			//////////////////////////////////////////////////////////////////////////////
			// Unlike the Overlap alignment, clamp the B cell values to zero
			tmp_b.sse = _MM_ADD(_MM_MAX(B_ptr->M.sse, zero.sse), gap_existance.sse);
			tmp_c.sse = _MM_ADD(_MM_MAX(B_ptr->I_target.sse, zero.sse), gap_extension.sse);
								
			X_ptr->I_target.sse = _MM_MAX(tmp_b.sse, tmp_c.sse);
			
			tmp_a.sse = _MM_CMPLT(tmp_b.sse, tmp_c.sse);
			
			X_ptr->I_target_start_i.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, B_ptr->M_start_i.sse),
							_mm_and_si128(
								tmp_a.sse, B_ptr->I_target_start_i.sse) );
			
			X_ptr->I_target_start_j.sse = _mm_or_si128(
							_mm_andnot_si128(
								tmp_a.sse, B_ptr->M_start_j.sse),
							_mm_and_si128(
								tmp_a.sse, B_ptr->I_target_start_j.sse) );
								
			//////////////////////////////////////////////////////////////////////////////
			// Unlike the Overlap alignment, find the maximum scoring element *anywhere*
			// in the dynamic programming matrix
			const sse_elem valid_query_and_target( 
				_mm_and_si128(valid_query.sse, _MM_CMPLT(all_j.sse, target_len.sse) ) );
			
			// tmp_a is true if this element is a new maximum M value
			tmp_a.sse = _mm_andnot_si128( _MM_CMPLT(X_ptr->M.sse, max_elem.M.sse), valid_query_and_target.sse );
			
			max_elem.M.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, X_ptr->M.sse),
						_mm_andnot_si128(tmp_a.sse, max_elem.M.sse) );
			
			max_elem.M_start_i.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, X_ptr->M_start_i.sse),
						_mm_andnot_si128(tmp_a.sse, max_elem.M_start_i.sse) );
			
			max_elem.M_start_j.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, X_ptr->M_start_j.sse),
						_mm_andnot_si128(tmp_a.sse, max_elem.M_start_j.sse) );
							
			stop_i.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, all_i.sse),
						_mm_andnot_si128(tmp_a.sse, stop_i.sse) );
			
			stop_j.sse = _mm_or_si128(
						_mm_and_si128(tmp_a.sse, all_j.sse),
						_mm_andnot_si128(tmp_a.sse, stop_j.sse) );			
		}
		
		// Swap the last_row and curr_row
		swap(last_row, curr_row);
	}
}

void SeqOverlap::init_BLOSUM62(vector<SO_Score> &m_matrix)
{
	// From the NCBI blast source distribution: ~/ncbi/data/BLOSUM62
	const SO_Score matrix[AA::MATRIX_SIZE] = {
		4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, -1,
		-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1,
		-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1,
		-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1,
		0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -1,
		-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1,
		-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1,
		0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1,
		-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1,
		-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1,
		-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1,
		-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1,
		-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1,
		-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1,
		-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -1,
		1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, -1,
		0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, -1,
		-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -1,
		-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1,
		0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1,
		-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1,
		-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1,
		-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
	};

	m_matrix = vector<SO_Score>(matrix, matrix + AA::MATRIX_SIZE);
}

base_type SO::na_to_bits(const char &m_base)
{
	switch(m_base){
		case 'A': case 'a':
			return NA::A;
		case 'C': case 'c':
			return NA::C;
		case 'G': case 'g':
			return NA::G;
		case 'T': case 't':
			return NA::T;
		case 'M': case 'm':
			return NA::M;
		case 'R': case 'r':
			return NA::R;
		case 'S': case 's':
			return NA::S;
		case 'V': case 'v':
			return NA::V;
		case 'W': case 'w': 
			return NA::W;
		case 'Y': case  'y':
			return NA::Y;
		case 'H': case 'h':
			return NA::H;
		case 'K': case 'k':
			return NA::K;
		case 'D': case 'd':
			return NA::D;
		case 'B': case 'b':
			return NA::B;
		case 'N': case 'n':
			return NA::N;
		case '-':
			return NA::GAP;
	};

	throw __FILE__ ":na_to_bits: Unknown base!";
	return NA::GAP; // Keep the compiler happy
};

base_type SO::aa_to_bits(const char &m_base)
{
	switch(m_base){
		case 'A': case 'a':
			return AA::A;
		case 'R': case 'r':
			return AA::R;
		case 'N': case 'n':
			return AA::N;
		case 'D': case 'd':
			return AA::D;
		case 'C': case 'c':
			return AA::C;
		case 'Q': case 'q':
			return AA::Q;
		case 'E': case 'e':
			return AA::E;
		case 'G': case 'g':
			return AA::G;
		case 'H': case 'h':
			return AA::H;
		case 'I': case 'i':
			return AA::I;
		case 'L': case 'l':
			return AA::L;
		case 'K': case 'k':
			return AA::K;
		case 'M': case 'm':
			return AA::M;
		case 'F': case 'f':
			return AA::F;
		case 'P': case 'p':
			return AA::P;
		case 'S': case 's':
			return AA::S;
		case 'T': case 't':
			return AA::T;
		case 'W': case 'w':
			return AA::W;
		case 'Y': case 'y':
			return AA::Y;
		case 'V': case 'v':
			return AA::V;
		case 'B': case 'b':
			return AA::B;
		case 'Z': case 'z':
			return AA::Z;
		case 'X': case 'x':
			return AA::X;
		case '-':
			return AA::GAP;
	};

	throw __FILE__ ":aa_to_bit: Unknown base!";
	return AA::GAP; // Keep the compiler happy
};
