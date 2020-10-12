#ifndef __JASONS_JSON
#define __JASONS_JSON

#include <unordered_map>
#include <deque>
#include <vector>
#include <string>
#include <algorithm>

#define		JSON_VERSION	"0.1"

namespace json
{
	class JSON
	{
		public:
			typedef enum {
				JSON_BOOL, 	// bool
				JSON_NUMBER, 	// double
				JSON_STRING, 	// std::string
				JSON_ARRAY, 	// std::deque
				JSON_MAP, 	// std::unordered_map
				JSON_NONE	// NULL
			} Type;
		
		// Only use a typedef for the more complicated
		// types (i.e. Array and Map). There is no need to replace std::string
		// with String (for now).
		typedef std::deque<JSON> Array;
		typedef std::unordered_map<std::string, JSON> Map;
		
		private:

			Type value_type;
			void* value_ptr;

			void assign(const std::string::const_iterator &m_begin, 
				const std::string::const_iterator &m_end);
			Type get_type(const std::string::const_iterator &m_begin, 
				const std::string::const_iterator &m_end);

			std::pair<void*, Type> parse(const std::string::const_iterator &m_begin, 
				const std::string::const_iterator &m_end);

			void* parse_none(const std::string::const_iterator &m_begin, 
				const std::string::const_iterator &m_end);
			void* parse_bool(const std::string::const_iterator &m_begin, 
				const std::string::const_iterator &m_end);
			void* parse_number(const std::string::const_iterator &m_begin, 
				const std::string::const_iterator &m_end);
			void* parse_string(const std::string::const_iterator &m_begin, 
				const std::string::const_iterator &m_end);
			void* parse_array(const std::string::const_iterator &m_begin, 
				const std::string::const_iterator &m_end);
			void* parse_map(const std::string::const_iterator &m_begin, 
				const std::string::const_iterator &m_end);

			double string_to_number(const std::string &m_str) const;
			
			void duplicate(const JSON& m_rhs)
			{
				clear();
				
				value_type = m_rhs.value_type;

				switch(value_type){
					case JSON_BOOL:

						value_ptr = new bool;

						if(value_ptr == NULL){
							throw __FILE__ ":JSON::duplicate: Unable to allocate bool";
						}

						*( (bool*) value_ptr ) = *( (bool*)m_rhs.value_ptr );
						break;
					case JSON_NUMBER:

						value_ptr = new double;

						if(value_ptr == NULL){
							throw __FILE__ ":JSON::duplicate: Unable to allocate double";
						}

						*( (double*) value_ptr ) = *( (double*) m_rhs.value_ptr );
						break;
					case JSON_STRING:

						value_ptr = new std::string;

						if(value_ptr == NULL){
							throw __FILE__ ":JSON::duplicate: Unable to allocate string";
						}

						*( (std::string*) value_ptr ) = *( (std::string*) m_rhs.value_ptr );

						break;
					case JSON_ARRAY:

						value_ptr = new std::deque<JSON>;

						if(value_ptr == NULL){
							throw __FILE__ ":JSON::duplicate: Unable to allocate array";
						}

						*( (Array*) value_ptr ) = *( (Array*)m_rhs.value_ptr );

						break;
					case JSON_MAP:

						value_ptr = new std::unordered_map<std::string, JSON>;

						if(value_ptr == NULL){
							throw __FILE__ ":JSON::duplicate: Unable to allocate map";
						}

						*( (Map*) value_ptr ) = 
							*( (Map*)m_rhs.value_ptr );

						break;
					case JSON_NONE:
						// Do nothing
						break;
					default:
						throw __FILE__ ":JSON::duplicate: Unknown type";
				};
			};
		public:
			
			JSON() : 
				value_type(JSON_NONE), value_ptr(NULL)
			{

			};

			JSON(const std::string &m_rhs) : 
				value_type(JSON_NONE), value_ptr(NULL)
			{
				assign( m_rhs.begin(), m_rhs.end() );
			};

			JSON(const JSON& m_rhs) : 
				value_type(JSON_NONE), value_ptr(NULL)
			{
				duplicate(m_rhs);
			};

			~JSON()
			{
				clear();
			};
			
			JSON& operator=(const std::string &m_rhs)
			{
				assign( m_rhs.begin(), m_rhs.end() );
				return *this;
			};

			JSON& operator=(const JSON &m_rhs)
			{
				JSON tmp;
				
				// Copy before clearing, since we may be copying from ourselves!
				tmp.duplicate(m_rhs);
				
				clear();
				
				// Move the memory in tmp to the current object
				value_type = tmp.value_type;
				value_ptr = tmp.value_ptr;
				
				// Make tmp safe to delete
				tmp.value_type = JSON_NONE;
				tmp.value_ptr = NULL;
				
				return *this;
			};

			void clear()
			{
				switch(value_type){
					case JSON_BOOL:

						if(value_ptr == NULL){
							throw __FILE__ ":JSON: Attempting to delete NULL bool ptr";
						};

						delete ( (bool*) value_ptr );
						break;
					case JSON_NUMBER:

						if(value_ptr == NULL){
							throw __FILE__ ":JSON: Attempting to delete NULL number ptr";
						};

						delete ( (double*) value_ptr );
						break;
					case JSON_STRING:

						if(value_ptr == NULL){
							throw __FILE__ ":JSON: Attempting to delete NULL string ptr";
						};

						delete ( (std::string*) value_ptr );

						break;
					case JSON_ARRAY:

						if(value_ptr == NULL){
							throw __FILE__ ":JSON: Attempting to delete NULL array ptr";
						};

						delete ( (Array*) value_ptr );

						break;
					case JSON_MAP:

						if(value_ptr == NULL){
							throw __FILE__ ":JSON: Attempting to delete NULL map ptr";
						};

						delete ( (Map*) value_ptr );

						break;
					case JSON_NONE:
						// Do nothing
						break;
					default:
						throw __FILE__ ":JSON: Unknown type";
				};

				value_ptr = NULL;
				value_type = JSON_NONE;
			};
			
			inline Type get_type() const
			{
				return value_type;
			};
			
			inline size_t size() const
			{
				if(value_ptr == NULL){
					throw __FILE__ ":JSON::size: NULL value_ptr";
				}
					
				if(value_type == JSON_MAP){
					return ( (Map*)value_ptr )->size();
				}
				
				if(value_type == JSON_ARRAY){
					return ( (Array*)value_ptr )->size();
				}
				
				throw __FILE__ ":JSON::size: Illegal call to non-map/non-array object";
				return 0;
			};
			
			inline bool get_bool(bool m_force_conversion = false) const
			{
				if(value_ptr == NULL){
					throw __FILE__ ":JSON::get_bool: NULL value_ptr";
				}
					
				if(value_type == JSON_BOOL){
					return *( (bool*)value_ptr );
				}
				
				// If the user wants to take a chance, we can convert number
				// to bool
				if(m_force_conversion){
					if(value_type == JSON_NUMBER){
						return ( *( (double*)value_ptr ) != 0 );
					}
				}
							
				throw __FILE__ ":JSON::get_bool: Non-bool value_ptr!";
				return false;
			};
			
			inline double get_number(bool m_force_conversion = false) const
			{
				if(value_ptr == NULL){
					throw __FILE__ ":JSON::get_number: NULL value_ptr";
				}
					
				if(value_type == JSON_NUMBER){					
					return *( (double*)value_ptr );
				}
				
				// If the user wants to take a chance, we can convert string and bool
				// to number
				if(m_force_conversion){
					
					if(value_type == JSON_BOOL){
						return double( *( (bool*)value_ptr ) );
					}
					
					if(value_type == JSON_STRING){
						return string_to_number( *( (std::string*)value_ptr ) );
					}
				}
				
				throw __FILE__ ":JSON::get_number: Non-number value_ptr!";
				return 0.0;
			};
			
			inline std::string get_string() const
			{
				if(value_type == JSON_STRING){
				
					if(value_ptr == NULL){
						throw __FILE__ ":JSON::get_string: NULL value_ptr";
					}
					
					return *( (std::string*)value_ptr );
				}
				
				throw __FILE__ ":JSON::get_string: Non-string value_ptr!";
				return std::string();
			};
			
			// Test for map element existance
			inline bool has_key(const std::string &m_key) const
			{
				if(value_type == JSON_MAP){
					
					Map::const_iterator iter = ( (Map*)value_ptr )->find(m_key);
					
					return ( iter != ( (Map*)value_ptr )->end() );
				}
				
				throw __FILE__ ":has_key: Attempted to dereference non-map object";
				return false; // Keep the compiler happy
			};
			
			inline std::vector<std::string> keys() const
			{
				if(value_ptr == NULL){
					throw __FILE__ ":JSON::keys: NULL value_ptr";
				}
					
				if(value_type == JSON_MAP){
					
					std::deque<std::string> local;
					
					for(Map::const_iterator i = ( (Map*)value_ptr )->begin();
						i != ( (Map*)value_ptr )->end();++i){
						
						local.push_back(i->first);
					}
					
					// Each Map key is unique, so there is no need to 
					// use std::unique on the extracted keys
					std::sort( local.begin(), local.end() );
					
					return std::vector<std::string>( local.begin(), local.end() );
				}
				
				throw __FILE__ ":JSON::keys: Illegal call to non-map object";
				return std::vector<std::string>();
			};
			
			// Access a map element
			const JSON& operator[](const std::string &m_key) const
			{
				if(value_type == JSON_MAP){
					
					Map::const_iterator iter = ( (Map*)value_ptr )->find(m_key);
					
					if( iter == ( (Map*)value_ptr )->end() ){
						throw __FILE__ ":JSON[map]: Did not find requested key in map";
					}
					
					return iter->second;
				}
				
				throw __FILE__ ":JSON[map]: Attempted to dereference non-map object";
				return *this; // Keep the compiler happy
			};
			
			// Access an array element
			const JSON& operator[](const size_t &m_index) const
			{
				if(value_type == JSON_ARRAY){
					
					if( m_index >= ( (Array*)value_ptr )->size() ){
						throw __FILE__ ":JSON[array]: Index out of bounds";
					}
					
					return (* ( (Array*)value_ptr ) )[m_index];
				}
				
				throw __FILE__ ":JSON[array]: Attempted to dereference non-array object";
				return *this; // Keep the compiler happy
			};
	};
}

#endif // __JASONS_JSON
