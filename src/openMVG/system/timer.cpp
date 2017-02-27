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

#include "openMVG/system/timer.hpp"


namespace openMVG {
namespace system {

Timer::Timer()
{
  reset();
}

void Timer::reset()
{
  start_ = std::chrono::high_resolution_clock::now();
}

double Timer::elapsed() const
{
  const auto end_ = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::seconds>(end_ - start_).count();
}

double Timer::elapsedMs() const
{
  const auto end_ = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::milliseconds>(end_ - start_).count();
}

std::ostream& operator << (std::ostream& str, const Timer& t)
{
  return str << t.elapsed() << " s elapsed";
}

} // namespace system
} // namespace openMVG
