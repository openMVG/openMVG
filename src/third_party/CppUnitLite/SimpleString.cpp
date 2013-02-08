

#include "SimpleString.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

SimpleString::SimpleString ()
: buffer(new char [1])
{
	buffer [0] = '\0';
}


SimpleString::SimpleString (const char *otherBuffer)
: buffer (new char [strlen (otherBuffer) + 1])
{
	strcpy (buffer, otherBuffer);
}

SimpleString::SimpleString (const SimpleString& other)
{
	buffer = new char [other.size() + 1];
	strcpy(buffer, other.buffer);
}


SimpleString & SimpleString::operator= (const SimpleString& other)
{
	delete buffer;
	buffer = new char [other.size() + 1];
	strcpy(buffer, other.buffer);	
	return *this;
}


char *SimpleString::asCharString () const
{
	return buffer;
}

int SimpleString::size() const
{
	return strlen (buffer);
}

SimpleString::~SimpleString ()
{
	delete [] buffer;
}


bool operator== (const SimpleString& left, const SimpleString& right)
{
	return !strcmp (left.asCharString (), right.asCharString ());
}
