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

#include <openMVG/system/timer.hpp>
#include <cmath>

#ifdef _WIN32
# include <windows.h>
#else
# include <sys/time.h>
#endif

namespace openMVG {
namespace system {

Timer::Timer()
{
#ifndef HAVE_CXX11_CHRONO
#ifdef _WIN32
  LARGE_INTEGER freq;
  if (!QueryPerformanceFrequency(&freq))
  {
    const char *msg = "Failed to initialize high resolution timer!";
    std::cerr << msg << std::endl;
    throw std::runtime_error(msg);
  }
  frequency_ = static_cast<double>(freq.QuadPart);
#endif
#endif
  reset();
}

void Timer::reset()
{
#ifdef HAVE_CXX11_CHRONO
  start_ = std::chrono::high_resolution_clock::now();
#else

#ifdef _WIN32
  LARGE_INTEGER li_start_;
  QueryPerformanceCounter(&li_start_);
  start_ = static_cast<double>(li_start_.QuadPart);
#else
  timeval start;
  gettimeofday(&start, nullptr);
  start_ = start.tv_sec + start.tv_usec * 1e-6;
#endif

#endif // HAVE_CXX11_CHRONO
}

double Timer::elapsed() const
{
#ifdef HAVE_CXX11_CHRONO
  return elapsedMs() / 1000.;
#else

  double elapsed_;
#ifdef _WIN32
  LARGE_INTEGER end_;
  QueryPerformanceCounter(&end_);
  elapsed_ = (static_cast<double>(end_.QuadPart) - start_) / frequency_;
#else
  timeval end;
  gettimeofday(&end, nullptr);
  const double end_ = end.tv_sec + end.tv_usec * 1e-6;
  elapsed_ = end_ - start_;
#endif
  return elapsed_;
#endif // HAVE_CXX11_CHRONO
}

double Timer::elapsedMs() const
{
#ifdef HAVE_CXX11_CHRONO
  const auto end_ = std::chrono::high_resolution_clock::now();
  return std::chrono::duration_cast<std::chrono::milliseconds>(end_ - start_).count();
#else
  return elapsed() * 1000.;
#endif // HAVE_CXX11_CHRONO
}

std::ostream& operator << (std::ostream& str, const Timer& t)
{
  return str << t.elapsed() << " s elapsed";
}

std::string prettyTime(double durationMs)
{
  std::string out;

  const auto msecs = fmod(durationMs, 1000);
  durationMs /= 1000.;
  const std::size_t secs = std::size_t(fmod(durationMs, 60));
  durationMs /= 60.;
  const std::size_t mins = std::size_t(fmod(durationMs, 60));
  durationMs /= 60.;
  const std::size_t hours = std::size_t(fmod(durationMs, 24));
  durationMs /= 24.;
  const std::size_t days = durationMs;

  bool printed_earlier = false;
  if(days >= 1)
  {
    printed_earlier = true;
    out += (std::to_string(days) + "d ");
  }
  if(printed_earlier || hours >= 1)
  {
    printed_earlier = true;
    out += (std::to_string(hours) + "h ");
  }
  if(printed_earlier || mins >= 1)
  {
    printed_earlier = true;
    out += (std::to_string(mins) + "m ");
  }
  if(printed_earlier || secs >= 1)
  {
    printed_earlier = true;
    out += (std::to_string(secs) + "s ");
  }
  if(printed_earlier || msecs >= 1)
  {
    printed_earlier = true;
    out += (std::to_string(msecs) + "ms");
  }
  return out;
}

} // namespace system
} // namespace openMVG
