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

namespace openMVG
{
/**
* @brief Namespace handling various functions and classes about system programming
*/
namespace system
{

/**
* @brief Timer class with microsecond accuracy.
*/
class Timer
{
  public:

    /**
    * @brief Default constructor
    */
    Timer();

    /**
    * @brief Reset the timer to zero.
    */
    void reset();

    /**
    * @brief Get the elapsed time in seconds.
    * @return Time elapsed in second
    */
    double elapsed() const;


    /**
    * @brief Get the elapsed time in milliseconds.
    * @return Time elapsed in milliseconds
    */
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

/**
* print the elapsed time
* @brief out Stream in which time is written
* @param tim Timer to output
* @return stream after write operation
*/
std::ostream& operator << ( std::ostream& out , const Timer& tim );

} // namespace system
} // namespace openMVG

#endif // OPENMVG_SYSTEM_TIMER_HPP

