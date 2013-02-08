///////////////////////////////////////////////////////////////////////////////
//
// TESTRESULT.H
// 
// A TestResult is a collection of the history of some test runs.  Right now
// it just collects failures.
// 
///////////////////////////////////////////////////////////////////////////////



#ifndef TESTRESULT_H
#define TESTRESULT_H

class Failure;

class TestResult
{
public:
					TestResult (); 
  virtual ~TestResult() {};
	virtual void	testsStarted ();
	virtual void	addFailure (const Failure& failure);
	virtual void	testsEnded ();
  
  int getFailureCount() const {return failureCount;}

private:
	int				failureCount;
};

#endif
