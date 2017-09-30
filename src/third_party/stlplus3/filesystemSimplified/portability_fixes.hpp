#ifndef STLPLUS_PORTABILITY_FIXES
#define STLPLUS_PORTABILITY_FIXES
////////////////////////////////////////////////////////////////////////////////

//   Author:    Andy Rushton
//   Copyright: (c) Southampton University 1999-2004
//              (c) Andy Rushton           2004 onwards
//   License:   BSD License, see ../docs/license.html

//   Contains work arounds for OS or Compiler specific problems to try to make
//   them look more alike

//   It is strongly recommended that this header be included as the first
//   #include in every source file

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Problem with MicroSoft defining two different macros to identify Windows
////////////////////////////////////////////////////////////////////////////////

#if defined(_WIN32) || defined(_WIN64) || defined(_WIN32_WCE)
#define MSWINDOWS
#endif

////////////////////////////////////////////////////////////////////////////////
// Problems with unnecessary or unfixable compiler warnings
////////////////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
// Microsoft Visual Studio
// shut up the following irritating warnings
//   4290 - C++ exception specification ignored
//   4996 - 'xxxx' was declared deprecated
#pragma warning(disable: 4290 4996)
#endif

#ifdef __BORLANDC__
// Borland
// Shut up the following irritating warnings
//   8026 - Functions with exception specifications are not expanded inline
//   8027 - Functions with xxx are not expanded inline
#pragma warn -8026
#pragma warn -8027
#endif

////////////////////////////////////////////////////////////////////////////////
// in some old system headers, std::max and std::min were unsafe - this is no longer the case
// I created extra template function definitions minimum/maximum that avoid all the problems above
// these are still supported for backwards compatibility

namespace stlplus
{
  template<typename T> const T& maximum(const T& l, const T& r) {return l > r ? l : r;}
  template<typename T> const T& minimum(const T& l, const T& r) {return l < r ? l : r;}
}

////////////////////////////////////////////////////////////////////////////////
// problems with missing functions
////////////////////////////////////////////////////////////////////////////////

#ifdef MSWINDOWS
unsigned sleep(unsigned seconds);
#else
#include <unistd.h>
#endif

////////////////////////////////////////////////////////////////////////////////
// Function for establishing endian-ness
////////////////////////////////////////////////////////////////////////////////
// Different machine architectures store data using different byte orders.
// This is referred to as Big- and Little-Endian Byte Ordering. 
//
// The issue is: where does a pointer to an integer type actually point?
//
// In both conventions, the address points to the left of the word but:
// Big-Endian - The most significant byte is on the left end of a word
// Little-Endian - The least significant byte is on the left end of a word
//
// Bytes are addressed left to right, so in big-endian order byte 0 is the
// msB, whereas in little-endian order byte 0 is the lsB. For example,
// Intel-based machines store data in little-endian byte order so byte 0 is
// the lsB.
//
// This function establishes byte order at run-time

namespace stlplus
{
  bool little_endian(void);
}

////////////////////////////////////////////////////////////////////////////////
#endif
