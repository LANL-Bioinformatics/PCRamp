#include "JSON.h"
#include <math.h>

// DEBUG
#include <iostream>

using namespace std;
using namespace json;

string::const_iterator trim_white_space_prefix(const string::const_iterator &m_begin, 
	const string::const_iterator& m_end);
string::const_iterator trim_white_space_suffix(const string::const_iterator &m_begin, 
	const string::const_iterator &m_end);
bool match_string(const string &m_str, 
	const string::const_iterator &m_begin, 
	const string::const_iterator &m_end);
string::const_iterator trim_delim_prefix(const char &m_delim, 
	const string::const_iterator &m_begin, 
	const string::const_iterator &m_end);
string::const_iterator trim_delim_suffix(const char &m_delim, 
	const string::const_iterator &m_begin, 
	const string::const_iterator &m_end);
string::const_iterator next(const char &m_delim, const string::const_iterator &m_begin, 
	const string::const_iterator &m_end);

void JSON::assign(const string::const_iterator &m_begin, 
				const string::const_iterator &m_end)
{
	clear();

	pair<void*, Type> ret = parse(m_begin, m_end);

	value_ptr = ret.first;
	value_type = ret.second;	
}

JSON::Type JSON::get_type(const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	string::const_iterator iter = trim_white_space_prefix(m_begin, m_end);

	if(iter == m_end){
		throw __FILE__ ":JSON::get_type: Unable to find type indicator";
	}

	switch(*iter){
		case '[':
			return JSON::JSON_ARRAY;
		case '{':
			return JSON::JSON_MAP;
		case '"':
			return JSON::JSON_STRING;
		case 'T': case 't': case 'F': case 'f':
			return JSON::JSON_BOOL;
		case '+': case '-': case '0': case '1': case '2': case '3':
		case '4': case '5': case '6': case '7': case '8': case '9':
			return JSON::JSON_NUMBER;
		case 'N': case 'n':
			return JSON::JSON_NONE;
		default:
			throw __FILE__ ":JSON::get_type: Unable to infer JSON type";
	}

	return JSON::JSON_NONE;
}

pair<void*, JSON::Type> JSON::parse(const string::const_iterator &m_begin, const string::const_iterator &m_end)
{
	pair<void*, JSON::Type> ret(NULL, JSON::JSON_NONE);

	ret.second = get_type(m_begin, m_end);

	switch(ret.second){
		case JSON::JSON_BOOL:
			ret.first = parse_bool(m_begin, m_end);
			break;
		case JSON::JSON_NUMBER:
			ret.first = parse_number(m_begin, m_end);
			break;
		case JSON::JSON_STRING:
			ret.first = parse_string(m_begin, m_end);
			break;
		case JSON::JSON_ARRAY:
			ret.first = parse_array(m_begin, m_end);
			break;
		case JSON::JSON_MAP:
			ret.first = parse_map(m_begin, m_end);
			break;
		case JSON::JSON_NONE:
			ret.first = parse_none(m_begin, m_end);
			break;
		default:
			throw __FILE__ ":JSON::parse: Unknown JSON type!";
			break;
	};

	return ret;
}

string::const_iterator trim_white_space_prefix(const string::const_iterator &m_begin, 
	const string::const_iterator& m_end)
{
	string::const_iterator ret = m_begin;

	while( (ret != m_end) && isspace(*ret) ){
		++ret;
	}

	if(ret == m_end){
		throw __FILE__ ":trim_white_space_prefix: Empty string";
	}

	return ret;
}

string::const_iterator trim_white_space_suffix(const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	if(m_begin == m_end){
		throw __FILE__ ":trim_white_space_suffix: Empty string";
	}

	string::const_iterator ret = m_end;

	while(ret != m_begin){

		string::const_iterator iter = ret - 1;

		if( isspace(*iter) ){
			ret = iter;
		}
		else{
			break;
		}
	}

	if(ret == m_begin){
		throw __FILE__ ":trim_white_space_suffix: Empty string";
	}

	return ret;
}

string::const_iterator trim_delim_prefix(const char &m_delim, 
	const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	string::const_iterator begin = trim_white_space_prefix(m_begin, m_end);
	
	if( (begin == m_end) || (*begin != m_delim) ){
		throw __FILE__ ":trim_delim_prefix: Unable to find delim";
	}

	// Skip the delimeter
	++begin;

	if(begin == m_end){
		throw __FILE__ ":trim_delim_prefix: Unable to find delim (2)";
	}

	return begin;
}

string::const_iterator trim_delim_suffix(const char &m_delim, 
	const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	string::const_iterator end = trim_white_space_suffix(m_begin, m_end);
	
	if( (m_begin == end) || ( *(end - 1) != m_delim ) ){
		throw __FILE__ ":trim_delim_suffix: Unable to find delim";
	}

	// Skip the delimeter
	--end;

	if(m_begin == end){
		throw __FILE__ ":trim_delim_suffix: Unable to find delim (2)";
	}

	return end;
}

bool match_string(const string &m_str, 
	const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	string::const_iterator iter = m_begin;

	for(string::const_iterator i = m_str.begin();i != m_str.end();++i, ++iter){

		if(iter == m_end){
			return false;
		}

		if( tolower(*i) != tolower(*iter) ){
			return false;
		}
	}

	return true;
}

string::const_iterator next(const char &m_delim, 
	const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	int square_count = 0;
	int curly_count = 0;
	int quote_count = 0;

	string::const_iterator ret = m_begin;

	while( (ret != m_end) ){

		if(*ret == m_delim){

			if( (square_count == 0) && (curly_count == 0) && 
				(quote_count%2 == 0) ){
				return ret;
			}
		}

		switch(*ret){
			case '[':
				++square_count;
				break;
			case ']':
				--square_count;
				break;
			case '{':
				++curly_count;
				break;
			case '}':
				--curly_count;
				break;
			case '"':
				++quote_count;
				break;
		};

		++ret;
	}

	return ret;
}

void* JSON::parse_none(const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	string::const_iterator begin = trim_white_space_prefix(m_begin, m_end);
	string::const_iterator end = trim_white_space_suffix(m_begin, m_end);

	if( !match_string("null", begin, end) ){
		throw __FILE__ ":: Unable to find null keyword";
	}

	return NULL;
};

void* JSON::parse_bool(const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	// DEBUG
	//cerr << "parse_bool" << endl;
	
	string::const_iterator begin = trim_white_space_prefix(m_begin, m_end);
	string::const_iterator end = trim_white_space_suffix(m_begin, m_end);

	bool *ret = new bool;

	if(ret == NULL){
		throw __FILE__ ":JSON::parse_bool: Unable to allocate boolean variable";
	}

	if( match_string("true", begin, end) ){

		*ret = true;
		return ret;
	}

	if( match_string("false", begin, end) ){

		*ret = false;
		return ret;
	}

	throw __FILE__ ":: Unable to find valid boolean keyword";
	return NULL;
};

void* JSON::parse_number(const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	// DEBUG
	//cerr << "parse_number" << endl;
	
	string::const_iterator begin = trim_white_space_prefix(m_begin, m_end);
	string::const_iterator end = trim_white_space_suffix(m_begin, m_end);

	double *ret = new double;

	if(ret == NULL){
		throw __FILE__ ":JSON::parse_number: Unable to allocate number";
	}

	const string local(begin, end);

	*ret = atof( local.c_str() );

	return ret;
}

void* JSON::parse_string(const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	// DEBUG
	//cerr << "parse_string" << endl;
	
	string::const_iterator begin = trim_delim_prefix('"', m_begin, m_end);
	string::const_iterator end = trim_delim_suffix('"', m_begin, m_end);

	string* ret = new string;

	if(ret == NULL){
		throw __FILE__ ":JSON::parse_string: Unable to allocate string";
	}

	ret->assign(begin, end);

	return ret;
}

void* JSON::parse_array(const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	// DEBUG
	//cerr << "parse_array" << endl;
	
	string::const_iterator begin = trim_delim_prefix('[', m_begin, m_end);
	string::const_iterator end = trim_delim_suffix(']', m_begin, m_end);

	Array* ret = new Array;

	if(ret == NULL){
		throw __FILE__ ":JSON::parse_array: Unable to allocate array";
	}

	while(true){

		string::const_iterator iter = next(',', begin, end);

		ret->push_back( JSON() );
		ret->back().assign(begin, iter);

		if(iter == end){
			break;
		}

		begin = iter + 1;
	}

	return ret;
}

void* JSON::parse_map(const string::const_iterator &m_begin, 
	const string::const_iterator &m_end)
{
	// DEBUG
	//cerr << "parse_map" << endl;
	
	string::const_iterator begin = trim_delim_prefix('{', m_begin, m_end);
	string::const_iterator end = trim_delim_suffix('}', m_begin, m_end);

	Map* ret = new Map;

	if(ret == NULL){
		throw __FILE__ ":JSON::parse_map: Unable to allocate array";
	}

	while(true){

		// DEBUG
		//cerr << "Looking for key" << endl;
		
		string::const_iterator iter = next(':', begin, end);

		if(iter == end){
			throw __FILE__ ":JSON::parse_map: Unable to find key/value delimeter";
		}

		string* key_ptr = (string*)parse_string(begin, iter);

		if(key_ptr == NULL){
			throw __FILE__ ":JSON::parse_map: Unable to parse key";
		}

		// DEBUG
		//cerr << "\tkey = " << *key_ptr << endl;
		
		begin = iter + 1;

		iter = next(',', begin, end);

		JSON value;

		// DEBUG
		//cerr << "value from: " << string(begin, iter) << endl;
		
		value.assign(begin, iter);

		// DEBUG
		//cerr << "About to copy" << endl;
		
		ret->insert( make_pair(*key_ptr, value) );
		
		delete key_ptr;
		key_ptr = NULL;

		if(iter == end){
			break;
		}

		begin = iter + 1;
	}
		
	return ret;
}

double JSON::string_to_number(const string &m_str) const
{
	// DEBUG
	//cerr << "string_to_number" << endl;
	
	string::const_iterator begin = trim_white_space_prefix( m_str.begin(), m_str.end() );
	string::const_iterator end = trim_white_space_suffix( m_str.begin(), m_str.end() );
	
	if(begin == end){
		throw __FILE__ ":JSON::string_to_number: Empty string";
	}
	
	double left_mantissa = 0.0; // Mantissa to the left of the decimal
	double right_mantissa = 0.0;// Mantissa to the right of the decimal
	double exponent = 0.0;
	double mantissa_sign = 1.0;
	double exponent_sign = 1.0;
	
	double decimal_pow = 1.0;
	
	enum{
		SIGN = 1,
		DECIMAL_POINT = 1 << 1,
		EXPONENT = 1 << 2,
		EXPONENT_SIGN = 1 << 3
	};
	
	unsigned int status = 0;
	
	for(;begin != end;++begin){
	
		switch(*begin){
			case '+':
				if(status & EXPONENT){
					if(status & EXPONENT_SIGN){
						throw __FILE__ ":JSON::string_to_number: Illegal exponent + sign";
					}
					else{
						status |= EXPONENT_SIGN;
					}
				}
				else{
					if(status & SIGN){
						throw __FILE__ ":JSON::string_to_number: Illegal + sign";
					}
					else{
						status |= SIGN;
					}
				}
				
				break;
			case '-':
				if(status & EXPONENT){
					if(status & EXPONENT_SIGN){
						throw __FILE__ ":JSON::string_to_number: Illegal exponent - sign";
					}
					else{
						status |= EXPONENT_SIGN;
						exponent_sign = -1.0;
					}
				}
				else{
					if(status & SIGN){
						throw __FILE__ ":JSON::string_to_number: Illegal - sign";
					}
					else{
						status |= SIGN;
						mantissa_sign = -1.0;
					}
				}
				
				break;
			case 'e': case 'E':
				if(status & EXPONENT){
					throw __FILE__ ":JSON::string_to_number: Multiple exponents";
				}
				else{
					status |= EXPONENT;
				}
				
				break;
			case '.':
				if(status & EXPONENT){
					throw __FILE__ ":JSON::string_to_number: Illegal . after exponent";
				}
				else{
					if(status & DECIMAL_POINT){
						throw __FILE__ ":JSON::string_to_number: Illegal multiple .";
					}
					
					status |= DECIMAL_POINT;
				}
				
				break;
			case '0': case '1': case '2': case '3': case '4':
			case '5': case '6': case '7': case '8': case '9':
				
				if(status & EXPONENT){
					exponent = 10.0*exponent + (*begin - '0');
				}
				else{ // Mantissa
					if(status & DECIMAL_POINT){ // Right of decimal
					
						decimal_pow /= 10.0;
						right_mantissa += decimal_pow*(*begin - '0');
					}
					else{ // Left of decimal
						left_mantissa = 10.0*left_mantissa + (*begin - '0');
					}
				}
				
				break;
			default:
				throw __FILE__ ":JSON::string_to_number: Illegal symbol";
			
		};
	}
	
	const double ret = mantissa_sign*(left_mantissa + right_mantissa)*pow(10.0, exponent_sign*exponent);
		
	return ret;
}
