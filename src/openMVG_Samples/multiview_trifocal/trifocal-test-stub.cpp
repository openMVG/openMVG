#include "testing/testing.h"


TEST(blabla, blele)
{
  CHECK(true);
}


/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */

