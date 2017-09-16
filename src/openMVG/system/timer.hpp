// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (C) 2013 David Ok <david.ok8@gmail.com>
// Copyright (C) 2014 Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_SYSTEM_TIMER_HPP
#define OPENMVG_SYSTEM_TIMER_HPP

#include <chrono>
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
    std::chrono::high_resolution_clock::time_point start_;
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
