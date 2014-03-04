////////////////////////////////////////////////////////////////////////////////

//   Author:    Andy Rushton
//   Copyright: (c) Southampton University 1999-2004
//              (c) Andy Rushton           2004 onwards
//   License:   BSD License, see ../docs/license.html

//   Notes:

//   Feb 2007: Rewritten in terms of platform-specific fixes to the
//   buffer-overflow problem. Using native functions for this has the added
//   benefit of giving access to other features of the C-runtime such as Unicode
//   support.

////////////////////////////////////////////////////////////////////////////////

#include "dprintf.hpp"
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <ctype.h>
#include <stdlib.h>

////////////////////////////////////////////////////////////////////////////////

namespace stlplus
{

////////////////////////////////////////////////////////////////////////////////

  int vdprintf(std::string& formatted, const char* format, va_list args)
  {
#ifdef MSWINDOWS
    int length = 0;
    char* buffer = 0;
    for(int buffer_length = 256; ; buffer_length*=2)
    {
      buffer = (char*)malloc(buffer_length);
      if (!buffer) return -1;
      length = _vsnprintf(buffer, buffer_length-1, format, args);
      if (length >= 0)
      {
        buffer[length] = 0;
        formatted += std::string(buffer);
        free(buffer);
        break;
      }
      free(buffer);
    }
    return length;
#else
    char* buffer = 0;
    int length = vasprintf(&buffer, format, args);
    if (!buffer) return -1;
    if (length >= 0)
      formatted += std::string(buffer);
    free(buffer);
    return length;
#endif
  }

  int dprintf(std::string& formatted, const char* format, ...)
  {
    va_list args;
    va_start(args, format);
    int result = vdprintf(formatted, format, args);
    va_end(args);
    return result;
  }

  std::string vdformat(const char* format, va_list args) throw(std::invalid_argument)
  {
    std::string formatted;
    int length = vdprintf(formatted, format, args);
    if (length < 0) throw std::invalid_argument("dprintf");
    return formatted;
  }

  std::string dformat(const char* format, ...) throw(std::invalid_argument)
  {
    std::string formatted;
    va_list args;
    va_start(args, format);
    int length = vdprintf(formatted, format, args);
    va_end(args);
    if (length < 0) throw std::invalid_argument("dprintf");
    return formatted;
  }

////////////////////////////////////////////////////////////////////////////////

} // end namespace stlplus
