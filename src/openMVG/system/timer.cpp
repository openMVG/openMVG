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
#ifdef _WIN32
# include <windows.h>
#else
# include <sys/time.h>
#endif

namespace openMVG {
namespace system {

  Timer::Timer()
   {
#ifdef HAVE_CXX11_CHRONO
#else
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
    gettimeofday(&start, NULL);
    start_ = start.tv_sec + start.tv_usec * 1e-6;
#endif

#endif // HAVE_CXX11_CHRONO
  }

  double Timer::elapsed() const
  {
    double elapsed_;
#ifdef HAVE_CXX11_CHRONO
    const auto end_ = std::chrono::high_resolution_clock::now();
    elapsed_ = std::chrono::duration_cast<std::chrono::milliseconds>(end_ - start_).count() / 1000.;
#else

#ifdef _WIN32
    LARGE_INTEGER end_;
    QueryPerformanceCounter(&end_);
    elapsed_ = (static_cast<double>(end_.QuadPart) - start_) / frequency_;
#else
    timeval end;
    gettimeofday(&end, NULL);
    const double end_ = end.tv_sec + end.tv_usec * 1e-6;
    elapsed_ = end_ - start_;
#endif
#endif // HAVE_CXX11_CHRONO
    return elapsed_;
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

} // namespace system
} // namespace openMVG
