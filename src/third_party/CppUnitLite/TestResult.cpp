
#include "TestResult.h"
#include "Failure.h"

#include <stdio.h>


TestResult::TestResult ()
	: failureCount (0)
{
}


void TestResult::testsStarted () 
{
}


void TestResult::addFailure (const Failure& failure) 
{
	fprintf (stdout, "%s%s%s%s%ld%s%s\n",
		"Failure: \"",
		failure.message.asCharString (),
		"\" " ,
		"line ",
		failure.lineNumber,
		" in ",
		failure.fileName.asCharString ());
		
	failureCount++;
}


void TestResult::testsEnded () 
{
  fprintf (stdout, "\n---------\n[TEST SESSION RESULTS]\n---------\n");
	if (failureCount > 0)
		fprintf (stdout, "There were %d failures\n", failureCount);
	else
		fprintf (stdout, "There were no test failures\n");
}
