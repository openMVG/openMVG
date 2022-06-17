// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (C) 2017 Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SYSTEM_PROGRESS_HPP
#define OPENMVG_SYSTEM_PROGRESS_HPP

#include <atomic>
#include <cstdint>
#include <string>

namespace openMVG {
namespace system {

//
// Class to keep track of the progress of a task.
// 1. Set the number of expected steps
// 2. Increment the number of completed step on the go
// 3. Ask for pourcentage along the process.
// Example:
// const std::uint32_t kExpectedCount = 300;
// ProgressInterface progress(kExpectedCount);
// for (std::uint32_t i = 0; i < kExpectedCount; ++i)
// {
//   std::cout << progress.Percent() << " ";
//   ++progress;
// }
//
class ProgressInterface {
 public:
  /** @brief Standard constructor
   * @param[in] expected_count The number of step of the process
   **/
  explicit ProgressInterface(const std::uint32_t expected_count = 1) noexcept
  {
    Restart(expected_count);
  }

  /// Destructor
  virtual ~ProgressInterface() = default;

  /** @brief Initializer of the ProgressInterface class
   * @param[in] expected_count The number of step of the process
   * @param[in] msg an optional status message
   * @return void if the progress class can be initialized, else it return false
   **/
  virtual void Restart(const std::uint32_t expected_count, const std::string& msg = {})
  {
    count_ = 0;
    expected_count_ = expected_count;
  }

  /** @brief Indicator if the current operation should be aborted.
   * @return Return true if the process has been canceled by the user.
   **/
  virtual bool hasBeenCanceled()const { return false; }

  /**
   * @brief Post-Increment operator
   * @param[in] increment the number of step that we want to increment the internal step counter
   * @return the value of the internal count => count()
   **/
  virtual std::uint32_t operator+=(const std::uint32_t increment)
  {
    count_ += increment;
    return count_;
  }

  /**
   * @brief Pre-Increment operator
   * @return the value of count_
   **/
  virtual std::uint32_t operator++()
  {
    return operator+=(1);
  }

  /** @brief A dummy progress reporter. Does nothing.
   **/
  static ProgressInterface& dummy()
  {
    static ProgressInterface null;
    return null;
  }

  /**
   * @brief Get the percent progress of the process
   * @return The percentage of the progress
   **/
  int Percent() const
  {
    return static_cast<int>(count_ / static_cast<float>(expected_count_) * 100.f + 0.5f);
  }

  //--
  //-- Accessor
  //--

  /**
   * @brief Get the internal count step
   * @return The current progress value
   **/
  std::uint32_t count() const { return count_; }

  /**
   * @brief Get the value of expected_count_
   * @return The value of expected_count_
   **/
  std::uint32_t expected_count() const { return expected_count_; }

 protected:
  /// Number of expected number of steps
  std::uint32_t expected_count_;
  /// Tracking of the number of completed steps => count_ will evolve in [0, expected_count_]
  std::atomic<std::uint32_t> count_;
};

} // namespace system
} // namespace openMVG
#endif // OPENMVG_SYSTEM_PROGRESS_HPP
