

#include "Test.h"
#include "TestRegistry.h"
#include "TestResult.h"
#include "Failure.h"


Test::Test (const SimpleString& testName) 
	: name_ (testName) 
{
	TestRegistry::addTest (this);
}



Test *Test::getNext() const
{
	return next_;
}

void Test::setNext(Test *test)
{	
	next_ = test;
}

bool Test::check(long expected, long actual, TestResult& result, const SimpleString& fileName, long lineNumber)
{
	if (expected == actual)
		return true;
	result.addFailure (
		Failure (
			name_, 
			StringFrom (__FILE__), 
			__LINE__, 
			StringFrom (expected), 
			StringFrom (actual)));

	return false;

}


bool Test::check(const SimpleString& expected, const SimpleString& actual, TestResult& result, const SimpleString& fileName, long lineNumber)
{
	if (expected == actual)
		return true;
	result.addFailure (
		Failure (
			name_, 
			StringFrom (__FILE__), 
			__LINE__, 
			expected, 
			actual));

	return false;

}

