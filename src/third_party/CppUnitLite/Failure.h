
///////////////////////////////////////////////////////////////////////////////
//
// FAILURE.H
//
// Failure is a class which holds information pertaining to a specific
// test failure.  The stream insertion operator is overloaded to allow easy
// display.
//
///////////////////////////////////////////////////////////////////////////////


#ifndef FAILURE_H
#define FAILURE_H

#include "SimpleString.h"


class Failure
{

public:
	Failure (const SimpleString&		theTestName, 
			 const SimpleString&		theFileName, 
			 long	  					theLineNumber,
			 const SimpleString&		theCondition);

	Failure (const SimpleString&		theTestName, 
			 const SimpleString&		theFileName, 
			 long						theLineNumber,
			 const SimpleString&		expected,
			 const SimpleString&		actual);

	SimpleString		message;
	SimpleString		testName;
	SimpleString		fileName;
	long				lineNumber;
};


#endif
