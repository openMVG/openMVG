#ifndef STLPLUS_TIMER
#define STLPLUS_TIMER
////////////////////////////////////////////////////////////////////////////////

//   Author:    Andy Rushton
//   Copyright: (c) Southampton University 1999-2004
//              (c) Andy Rushton           2004 onwards
//   License:   BSD License, see ../docs/license.html

//   A CPU timer encapsulated as a class. Measures the CPU time used since its
//   construction and allows this cumulative time to be reported at any time.

////////////////////////////////////////////////////////////////////////////////
#include "subsystems_fixes.hpp"
#include <time.h>
#include <string>
#include <iostream>

namespace stlplus
{

  ////////////////////////////////////////////////////////////////////////////////

  class timer
  {
  private:
    clock_t m_clock;
    time_t m_time;

  public:
    // constructor resets the timer to zero
    timer(void);
    ~timer(void);

    // reset the timer to zero without destroying it
    void reset(void);

    // get the elapsed time in seconds, expressed as a float
    float elapsed(void) const;
    // get the CPU time in seconds, expressed as a float
    float cpu(void) const;

    // get a printable string representing the elapsed time and CPU time
    //std::string text(void) const;
  };

  ////////////////////////////////////////////////////////////////////////////////

  // print the elapsed time and CPU time using the same representation as the text method
  std::ostream& operator << (std::ostream&, const timer&);

  ////////////////////////////////////////////////////////////////////////////////

} // end namespace stlplus

#endif
