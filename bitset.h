#ifndef __BITSET
#define __BITSET

#include <vector>

// Extend std::vector<bool> with set operations and MPI serialization
class BitSet : public std::vector<bool>
{
	public:
		
		BitSet()
		{
		};
		
		BitSet(const unsigned int &m_len, const bool &m_value)
		{
			resize(m_len, m_value);
		};
		
		inline unsigned int count() const
		{
			unsigned int ret = 0;
			
			for(const_iterator i = begin();i != end();++i){
				
				if(*i){
					++ret;
				}
			}
			
			return ret;
		};
		
		// Set union
		inline BitSet operator|(const BitSet &m_rhs) const
		{
			const size_t len = size();
			
			assert( len == m_rhs.size() );
			
			BitSet ret(len, false);
			
			for(size_t i = 0;i < len;++i){
				ret[i] = ( (*this)[i] || m_rhs[i] );
			}
			
			return ret;
		};
		
		// Set union
		inline BitSet& operator|=(const BitSet &m_rhs)
		{
			const size_t len = size();
			
			assert( len == m_rhs.size() );			
			
			for(size_t i = 0;i < len;++i){
				(*this)[i] = ( (*this)[i] || m_rhs[i] );
			}
			
			return *this;
		};
		
		// Set intersection
		inline BitSet operator&(const BitSet &m_rhs) const
		{
			const size_t len = size();
			
			assert( len == m_rhs.size() );
			
			BitSet ret(len, false);
			
			for(size_t i = 0;i < len;++i){
				ret[i] = ( (*this)[i] && m_rhs[i] );
			}
			
			return ret;
		};
		
		// Set intersection
		inline BitSet& operator&=(const BitSet &m_rhs)
		{
			const size_t len = size();
			
			assert( len == m_rhs.size() );			
			
			for(size_t i = 0;i < len;++i){
				(*this)[i] = ( (*this)[i] && m_rhs[i] );
			}
			
			return *this;
		};
		
		// Asymmetric set difference
		inline BitSet operator-(const BitSet &m_rhs) const
		{
			const size_t len = size();
			
			assert( len == m_rhs.size() );
			
			BitSet ret(*this);
			
			for(size_t i = 0;i < len;++i){
				
				if(m_rhs[i]){
					ret[i] = false;
				}
			}
			
			return ret;
		};
		
		inline bool empty_set() const
		{
			for(const_iterator i = begin();i != end();++i){
				
				if(*i){
					return false;
				}
			}
			
			return true;
		};
		
		template<class T> friend size_t mpi_size(const T &m_obj);
		template<class T> friend unsigned char* mpi_unpack(unsigned char* m_ptr, 
			const T &m_obj);
		template<class T> friend unsigned char* mpi_pack(unsigned char* m_ptr,
			const T &m_obj);
};

template<> size_t mpi_size(const BitSet &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const BitSet &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, BitSet &m_obj);

#endif // __BITSET
