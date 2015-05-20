// ========================================================================== //
//
// Copyright (C) 2013 David Ok <david.ok8@gmail.com>
// Copyright (C) 2014 Pierre Moulon
//
// Adapted from DO++, a basic set of libraries in C++ for computer 
// vision.
//
// This Source Code Form is subject to the terms of the Mozilla Public 
// License v. 2.0. If a copy of the MPL was not distributed with this file, 
// you can obtain one at http://mozilla.org/MPL/2.0/.
// ========================================================================== //

#ifndef OPENMVG_SYSTEM_TIMER_HPP
#define OPENMVG_SYSTEM_TIMER_HPP

#ifdef HAVE_CXX11_CHRONO
#include <chrono>
#endif
#include <iostream>

namespace openMVG {
namespace system {

  //! \brief Timer class with microsecond accuracy.
  class Timer
  {
  public:
    //! Default constructor
    Timer();
    //! Reset the timer to zero.
    void reset();
    //! Returns the elapsed time in seconds.
    double elapsed() const;
    //! Returns the elapsed time in milliseconds.
    double elapsedMs() const;
  private:

#ifdef HAVE_CXX11_CHRONO
    std::chrono::high_resolution_clock::time_point start_;
#else
    double start_;
#ifdef _WIN32
    double frequency_;
#endif
#endif // HAVE_CXX11_CHRONO
  };
  
  // print the elapsed time
  std::ostream& operator << (std::ostream&, const Timer&);

} // namespace system
} // namespace openMVG

#endif // OPENMVG_SYSTEM_TIMER_HPP

