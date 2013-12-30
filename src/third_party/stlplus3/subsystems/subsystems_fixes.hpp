#ifndef STLPLUS_SUBSYSTEMS_FIXES
#define STLPLUS_SUBSYSTEMS_FIXES
////////////////////////////////////////////////////////////////////////////////

//   Author:    Andy Rushton
//   Copyright: (c) Southampton University 1999-2004
//              (c) Andy Rushton           2004 onwards
//   License:   BSD License, see ../docs/license.html

//   Contains work arounds for OS or Compiler specific problems

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Unnecessary compiler warnings
////////////////////////////////////////////////////////////////////////////////

#ifdef _MSC_VER
// Microsoft Visual Studio
// shut up the following irritating warnings
//   4786 - VC6, identifier string exceeded maximum allowable length and was truncated (only affects debugger)
//   4305 - VC6, identifier type was converted to a smaller type
//   4503 - VC6, decorated name was longer than the maximum the compiler allows (only affects debugger)
//   4309 - VC6, type conversion operation caused a constant to exceeded the space allocated for it
//   4290 - VC6, C++ exception specification ignored
//   4800 - VC6, forcing value to bool 'true' or 'false' (performance warning)
//   4675 - VC7.1, "change" in function overload resolution _might_ have altered program
//   4996 - VC8, 'xxxx' was declared deprecated
#pragma warning(disable: 4786 4305 4503 4309 4290 4800 4675 4996)
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
#endif
