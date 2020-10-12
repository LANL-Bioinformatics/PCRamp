#ifndef __READ_ONLY_MULTIMAP
#define __READ_ONLY_MULTIMAP

#include <deque>
#include <algorithm>

#include "sort.h"

// A very simple, "read-only" multimap implementation that:
//	a) Saves space by storing the {key, value} pairs in a std::deque
//	b) Saves time by requiring an explicit (but parallel) call to the
//	   read_only_multimap::sort() member function to enable access (via
//	   begin and/or equal_range). The sorting member function must
//	   be called some time after the last insertion and before the first
//	   access.
//
// <no free lunch> To guarantee data integrity and ordering, this data structure
// does not allow modification of existing {key, value} pairs. </no free lunch>

template<class A, class B>
class read_only_multimap
{
	private:
		std::deque< std::pair<A, B> > data;
		bool ordered;
		
		struct key_comp
		{
			inline bool operator()( // for lower_bound
				const std::pair<A, B> &m_value,
				const A &m_key) const
			{
				return m_value.first < m_key;
			}
			
			inline bool operator()( // for upper_bound
				const A &m_key,
				const std::pair<A, B> &m_value) const
			{
				return m_key < m_value.first;
			}
			
			inline bool operator()( // for sort
				const std::pair<A, B> &m_lhs,
				const std::pair<A, B> &m_rhs) const
			{
				return m_lhs.first < m_rhs.first;
			}
		};

	public:

		typedef typename std::deque< std::pair<A, B> >::const_iterator const_iterator;
		
		read_only_multimap()
		{
			ordered = true;
		};
		
		inline void insert(const std::pair<A, B> &m_elem)
		{
			ordered = false;
			data.push_back(m_elem);
		};
		
		inline void insert(const const_iterator &m_begin, const const_iterator &m_end)
		{
			ordered = false;
			
			for(const_iterator i = m_begin;i != m_end;++i){
				data.push_back(*i);
			}
		};
		
		inline const_iterator begin() const
		{
			if(!ordered){
				
				// Since the begin() function is not allowed to make any modifications,
				// all we can do is throw an error (to let the programmer know that they
				// need to envoke the sort() member function before calling begin().
				throw __FILE__ ":read_only_multimap::begin: Cannot iterate an unsorted multimap!";
			}
			
			return data.begin();
		};
		
		inline const_iterator end() const
		{
			return data.end();
		};
		
		void sort()
		{
			if(!ordered){
				
				SORT( data.begin(), data.end(), key_comp() );
				ordered = true;
			}
		};
		
		std::pair<const_iterator, const_iterator> equal_range(const A &m_key) const
		{
			if(!ordered){
				
				// Since the begin() function is not allowed to make any modifications,
				// all we can do is throw an error (to let the programmer know that they
				// need to envoke the sort() member function before calling equal_range().
				throw __FILE__ ":read_only_multimap::equal_range: Cannot iterate an unsorted multimap!";
			}
			
			const_iterator iter = std::lower_bound( data.begin(), data.end(), m_key, key_comp() );
			
			return std::make_pair(iter, 
				std::upper_bound(iter, data.end(), m_key, key_comp() ) );
		};
		
		inline bool empty() const
		{
			return data.empty();
		};
		
		inline size_t size() const
		{
			return data.size();
		};
};

#endif // __READ_ONLY_MULTIMAP
