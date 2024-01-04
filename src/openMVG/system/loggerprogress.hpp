// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (C) 2017 Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SYSTEM_LOGPROGRESS_HPP
#define OPENMVG_SYSTEM_LOGPROGRESS_HPP

#include <openMVG/system/progressinterface.hpp>
#include <openMVG/system/logger.hpp>

#include <mutex>

namespace openMVG {
namespace system {

// Typical usage
// LoggerProgress progress(100, "My Task name", kModuloDisplay);
// for (int i = 0; i < 100; ++i) {
//   ++progress; // Will log directly to the LOG interface something like [My Task name] 10%
// }
//
// If you want to have a fine control over the LOG you can use
// progress.increment(1) that return a string. The string is empty if there is not need to display progress update:
//  (progress status modulo condition not meet)
// for (int i = 0; i < 100; ++i) {
//   const bool is_display_required = progress.Increment(1);
//   LOG_IF(INFO, is_display_required) << progress.PercentString(); // here the logger will write the actual code line number and file
// }

class LoggerProgress : public ProgressInterface {
 public:
  /** @brief Constructor
   * @param[in] expected_count The number of step of the process
   * @param[in] msg an optional status message that acts as the task name
   * @param[in] modulo_percentage Indicate how often you want to display the progress status. Value must be in ]0;100[
   **/
  LoggerProgress(
      const std::uint32_t expected_count=1,
      const std::string& msg = {},
      const int modulo_percentage = 10) noexcept
      : ProgressInterface(expected_count),
        msg_(""),
        previously_displayed_percentage_(0),
        modulo_percentage_(modulo_percentage)
  {
    if (modulo_percentage < 1 || modulo_percentage >= 100)
      OPENMVG_LOG_ERROR << "modulo_percentage must be within the range ]0;100[.";

    Restart(expected_count, msg);
  }

  /** @brief Initializer of the LoggerProgress class
   * @param[in] expected_count The number of step of the process
   * @param[in] msg an optional status message
   * @return true if the progress class can be initialized, else it return false
   **/
  void Restart(const std::uint32_t expected_count, const std::string& msg = {}) override
  {
    ProgressInterface::Restart(expected_count);
    msg_ = msg.empty() ? "" : "[" + msg + "]";
  }

  /**
   * @brief Post-Increment operator (Display the progress status only if needed)
   * @param[in] increment the number of step that we want to increment the internal step counter
   * @return the value of the internal count => count()
   **/
  inline std::uint32_t operator+=(const std::uint32_t increment) override
  {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto res = ProgressInterface::operator+=(increment);
    const auto percentage = Percent();
    // Check if we need to display the progress status or not
    if (percentage % modulo_percentage_ == 0 && previously_displayed_percentage_ != percentage) {
      OPENMVG_LOG_INFO << PercentString();
      previously_displayed_percentage_ = percentage;
    }
    return res;
  }

  /**
   * @brief Increment operator (Return the progress status only if needed)
   * @param[in] increment the number of step that we want to increment the internal step counter
   * @return true if the progress status must be displayed. Else it is not necessary to display it
   **/
  inline bool Increment(const int increment = 1)
  {
    std::ostringstream os;
    std::lock_guard<std::mutex> lock(mutex_);
    ProgressInterface::operator+=(increment);
    const auto percentage = Percent();
    // Check if we need to display this progress status or not
    if (percentage % modulo_percentage_ == 0 && previously_displayed_percentage_ != percentage) {
      previously_displayed_percentage_ = percentage;
      return true;
    }
    return false;
  }

  /**
   * @brief Pre-Increment operator (Display the progress status only if needed)
   * @return the value of the internal count => count()
   **/
   inline std::uint32_t operator++() override
   {
     return this->operator+=(1);
   }

  /// Return the percentage of process as a string of type "[msg_] X%"
  inline std::string PercentString() const
  {
    return msg_ + " " + std::to_string(Percent()) + "%";
  }

 private:
  /// Store the task name
  std::string msg_;

  /// Store the previously displayed value: used in order to avoid the same percentage displayed many time
  int previously_displayed_percentage_;
  /// Store how many percent we need to display the progress status
  int modulo_percentage_;
  /// Mutex used to sync status update of the various internal variable
  std::mutex mutex_;
};

} // namespace system
} // namespace openMVG
#endif // OPENMVG_SYSTEM_LOGPROGRESS_HPP
