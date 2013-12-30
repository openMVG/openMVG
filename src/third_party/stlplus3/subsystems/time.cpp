////////////////////////////////////////////////////////////////////////////////

//   Author:    Andy Rushton
//   Copyright: (c) Southampton University 1999-2004
//              (c) Andy Rushton           2004 onwards
//   License:   BSD License, see ../docs/license.html

////////////////////////////////////////////////////////////////////////////////
#include "time.hpp"
#include "dprintf.hpp"
#include <ctype.h>
////////////////////////////////////////////////////////////////////////////////

namespace stlplus
{

  ////////////////////////////////////////////////////////////////////////////////

  time_t time_now(void)
  {
    return time(0);
  }

  time_t localtime_create(int year, int month, int day, int hour, int minute, int second)
  {
    tm tm_time;
    tm_time.tm_year = year-1900;  // years are represented as an offset from 1900, for reasons unknown
    tm_time.tm_mon = month-1;     // internal format represents month as 0-11, but it is friendlier to take an input 1-12
    tm_time.tm_mday = day;
    tm_time.tm_hour = hour;
    tm_time.tm_min = minute;
    tm_time.tm_sec = second;
    tm_time.tm_isdst = -1;        // specify that the function should work out daylight savings
    time_t result = mktime(&tm_time);
    return result;
  }

  int localtime_year(time_t t)
  {
    tm* tm_time = localtime(&t);
    return tm_time->tm_year + 1900;
  }

  int localtime_month(time_t t)
  {
    tm* tm_time = localtime(&t);
    return tm_time->tm_mon + 1;
  }

  int localtime_day(time_t t)
  {
    tm* tm_time = localtime(&t);
    return tm_time->tm_mday;
  }

  int localtime_hour(time_t t)
  {
    tm* tm_time = localtime(&t);
    return tm_time->tm_hour;
  }

  int localtime_minute(time_t t)
  {
    tm* tm_time = localtime(&t);
    return tm_time->tm_min;
  }

  int localtime_second(time_t t)
  {
    tm* tm_time = localtime(&t);
    return tm_time->tm_sec;
  }

  int localtime_weekday(time_t t)
  {
    tm* tm_time = localtime(&t);
    return tm_time->tm_wday;
  }

  int localtime_yearday(time_t t)
  {
    tm* tm_time = localtime(&t);
    return tm_time->tm_yday;
  }

  std::string localtime_string(time_t t)
  {
    tm* local = localtime(&t);
    std::string result = local ? asctime(local) : "*time not available*";
    // ctime appends a newline for no apparent reason - clean up
    while (!result.empty() && isspace(result[result.size()-1]))
      result.erase(result.size()-1,1);
    return result;
  }

  std::string delaytime_string(time_t seconds)
  {
    unsigned minutes = (unsigned)seconds / 60;
    seconds %= 60;
    unsigned hours = minutes / 60;
    minutes %= 60;
    unsigned days = hours / 24;
    hours %= 24;
    unsigned weeks = days / 7;
    days %= 7;
    std::string result;
    if (weeks > 0)
      result += dformat("%dw ",weeks);
    if (!result.empty() || days > 0)
      result += dformat("%dd ", days);
    if (!result.empty() || hours > 0)
      result += dformat("%d:", hours);
    if (!result.empty() || minutes > 0)
    {
      if (!result.empty())
        result += dformat("%02d:", minutes);
      else
        result += dformat("%d:", minutes);
    }
    if (!result.empty())
      result += dformat("%02d:", seconds);
    else
      result += dformat("%ds", seconds);
    return result;
  }

  ////////////////////////////////////////////////////////////////////////////////

} // end namespace stlplus
