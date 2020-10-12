#ifndef __UPDATE
#define __UPDATE

#include <string>
#include <sstream>
#include <iostream>

class UpdateInfo : public std::stringstream {
	
	private:
		size_t buffer_size;
		std::ostream &out;
	public:
	
		UpdateInfo(const std::string &m_prefix, std::ostream &m_out = std::cerr) :
			out(m_out)
		{
			init(m_prefix);
		};
		
		UpdateInfo(std::ostream &m_out = std::cerr) : out(m_out)
		{
			init( std::string() );
		};
		
		void init(const std::string &m_prefix);
		void flush();
		void close();
};

#endif // __UPDATE
