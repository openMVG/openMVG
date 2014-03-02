#ifndef STLPLUS_TIME
#define STLPLUS_TIME
////////////////////////////////////////////////////////////////////////////////

//   Author:    Andy Rushton
//   Copyright: (c) Southampton University 1999-2004
//              (c) Andy Rushton           2004 onwards
//   License:   BSD License, see ../docs/license.html

//   Simplified access to representations of time and conversions between them.
//   The motivation for this package is that the low-level system calls for
//   accessing time are ugly and therefore potentially error-prone. I hope that
//   this interface is much simpler and therefore easier to use and more likely
//   to yield first-time right programs.

//   time is represented as the built-in integer type time_t - this is the
//   standard representation of system time in computerland and represents the
//   number of seconds since midnight 1 Jan 1970, believe it or not.

//   Functions are provided here for converting to and from more
//   human-comprehendable forms.

////////////////////////////////////////////////////////////////////////////////
#include "portability_fixes.hpp"
#include <string>
#include <time.h>

namespace stlplus
{

  // get the integer representing the time now
  time_t time_now(void);

  // get the integer representing the requested time - the local time is expressed in the local timezone
  time_t localtime_create(int year, int month, int day, int hour, int minute, int second);

  // extract human-centric form of the machine representation time_t
  int localtime_year(time_t);    // the year e.g. 1962
  int localtime_month(time_t);   // the month, numbered 1-12 e.g. August = 8
  int localtime_day(time_t);     // the day of the month numbered 1-31 e.g. 29
  int localtime_hour(time_t);    // the hour of day numbered 0-23
  int localtime_minute(time_t);  // minute past the hour numbered 0-59
  int localtime_second(time_t);  // second past the minute numbered 0-59
  int localtime_weekday(time_t); // the day of the week numbered 0-6 with 0=Sunday
  int localtime_yearday(time_t); // the number of days into the year

  // convert the integer representation of time to a human-readable form
  std::string localtime_string(time_t);

  // convert a time delay in seconds to human-readable form
  std::string delaytime_string(time_t);

} // end namespace stlplus

#endif
