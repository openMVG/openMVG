////////////////////////////////////////////////////////////////////////////////

//   Author:    Andy Rushton
//   Copyright: (c) Southampton University 1999-2004
//              (c) Andy Rushton           2004 onwards
//   License:   BSD License, see ../docs/license.html

////////////////////////////////////////////////////////////////////////////////
#include "timer.hpp"
#include "dprintf.hpp"
#include "time.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace stlplus
{

////////////////////////////////////////////////////////////////////////////////

  timer::timer(void)
  {
    reset();
  }

  timer::~timer(void)
  {
  }

  void timer::reset(void)
  {
    m_clock = clock();
    m_time = time(0);
  }

  float timer::cpu(void) const
  {
    return ((float)(clock() - m_clock)) / ((float)CLOCKS_PER_SEC);
  }

  float timer::elapsed(void) const
  {
    return difftime(time(0), m_time);
  }

  /*std::string timer::text(void) const
  {
    return dformat("%4.2fs CPU, %s elapsed", cpu(), delaytime_string(time(0)-m_time).c_str());
  }*/

////////////////////////////////////////////////////////////////////////////////

  std::ostream& operator << (std::ostream& str, const timer& t)
  {
    //return str << t.text();
    return str << t.cpu() << " s elapsed";
  }

////////////////////////////////////////////////////////////////////////////////

} // end namespace stlplus
