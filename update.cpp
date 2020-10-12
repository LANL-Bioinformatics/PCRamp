#include "update.h"
#include <iostream>

using namespace std;

void UpdateInfo::init(const string &m_prefix)
{
	out << m_prefix;
	
	buffer_size = 0;
}

void UpdateInfo::flush()
{
	// Clear the buffer
	for(size_t i = 0;i < buffer_size;i++){
		out << '\b';
	}
	
	for(size_t i = 0;i < buffer_size;i++){
		out << ' ';
	}
	
	for(size_t i = 0;i < buffer_size;i++){
		out << '\b';
	}
	
	const string tmp = str();
	
	buffer_size = tmp.size();
	
	out << tmp;
	
	// Clear the stringstream buffer
	str( string() );
}

void UpdateInfo::close()
{
	out << endl;
	
	// Clear the stringstream buffer
	str( string() );
}
