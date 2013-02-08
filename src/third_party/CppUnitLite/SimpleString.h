
///////////////////////////////////////////////////////////////////////////////
//
// SIMPLESTRING.H
//
// One of the design goals of CppUnitLite is to compilation with very old C++
// compilers.  For that reason, I've added a simple string class that provides
// only the operations needed in CppUnitLite.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef SIMPLE_STRING
#define SIMPLE_STRING

#include <iostream>
#include <sstream>
#include <stdio.h>

class SimpleString
{
	friend bool	operator== (const SimpleString& left, const SimpleString& right);

public:
						SimpleString ();
						SimpleString (const char *value);
						SimpleString (const SimpleString& other);
						~SimpleString ();

	SimpleString	&	operator= (const SimpleString& other);

	char				*asCharString () const;
	int					size() const;

private:
	char				*buffer;
};

template <class T>
SimpleString StringFrom (const T & other)
{
  std::ostringstream os;
  os << other;
  return SimpleString(os.str().c_str());
}
template <>
inline SimpleString StringFrom (const bool & value)
{
  std::ostringstream os;
  os << (value ? "true" : "false");
	return SimpleString(os.str().c_str());
}


#endif
