

#include "Test.h"
#include "TestResult.h"
#include "TestRegistry.h"
#include <typeinfo>

void TestRegistry::addTest (Test *test) 
{
	instance ().add (test);
}


int TestRegistry::runAllTests (TestResult& result) 
{
	return instance ().run (result);
}


TestRegistry& TestRegistry::instance () 
{
	static TestRegistry registry;
	return registry;
}


void TestRegistry::add (Test *test) 
{
	if (tests == 0) {
		tests = test;
		return;
	}
	
	test->setNext (tests);
	tests = test;
}

int TestRegistry::run (TestResult& result) 
{
	result.testsStarted ();
  int failureCountCounter = 0;
	for (Test *test = tests; test != 0; test = test->getNext ())  {
    std::cout<< "\n\n Launching Test : "<< typeid(*test).name()
      << std::endl << std::endl;
    test->run (result);
    if (failureCountCounter != result.getFailureCount())
    {
      failureCountCounter = result.getFailureCount();
      std::cout << "\n\n \t TEST -> [FAILURE]\n" << std::endl;
    }
    else
    {
      std::cout << "\n\n \t TEST -> [OK]\n" << std::endl;
    }
  }
	result.testsEnded ();
  return result.getFailureCount();
}



