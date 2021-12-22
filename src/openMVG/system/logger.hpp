// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (C) 2017 Pierre Moulon

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef OPENMVG_SYSTEM_LOGGER_HPP
#define OPENMVG_SYSTEM_LOGGER_HPP

#include <iostream>
#include <sstream>
#include <string>

namespace openMVG {
namespace system {
namespace logger {

enum class ELogMode : unsigned char
{
  VERBOSITY_INFO = 0,
  VERBOSITY_WARNING,
  VERBOSITY_ERROR
};

// Configure a default severity for the logger
static ELogMode logger_severity = ELogMode::VERBOSITY_INFO;

inline std::string ELogModeToString(const ELogMode & log_mode)
{
  switch (log_mode)
  {
    case ELogMode::VERBOSITY_INFO:
      return {"INFO: "};
    case ELogMode::VERBOSITY_WARNING:
      return {"WARNING: "};
    case ELogMode::VERBOSITY_ERROR:
      return {"ERROR: "};
    default:
      return {"UNKNOWN: "};
  }
  return {"UNKNOWN: "};
}

class StreamMessageLogger
{
private:
  std::ostringstream ostream_;

public:
  StreamMessageLogger() = default;

  //--
  // Make this class non copyable
  StreamMessageLogger( const StreamMessageLogger & other ) = delete; // non construction-copyable
  StreamMessageLogger & operator=( const StreamMessageLogger & ) = delete; // non copyable

  ~StreamMessageLogger()
  {
    ostream_ << '\n';
    // Output the contents of the stream (cerr by default)
    std::cerr << ostream_.str();
  }

  std::ostringstream & ostream() { return ostream_; }

  StreamMessageLogger & setMessage(
    ELogMode mode = ELogMode::VERBOSITY_INFO,
    const std::string & file = "",
    const std::string & line = "",
    const std::string & message = "")
  {
    ostream_
      << ELogModeToString(mode)
      << '['<< file << ':' << line << "] "
      << message;
    return *this;
  }
};

class NullBuffer : public std::streambuf
{
public:
  void operator&(const std::ostream &s) { }
};

// strips filepath to get filename (part after / or \ if any)
inline const char* filename(const char* path)
{
  for (auto ptr = path; *ptr; ++ptr)
  {
    if (*ptr == '/' || *ptr == '\\')
    {
      path = ptr + 1;
    }
  }
  return path;
}

} // namespace logger
} // namespace system
} // namespace openMVG

#define OPENMVG_LOG(MODE) \
  (MODE < openMVG::system::logger::logger_severity)\
    ? (void) 0\
    : openMVG::system::logger::NullBuffer() &\
        openMVG::system::logger::StreamMessageLogger().setMessage(\
          MODE,\
          openMVG::system::logger::filename(__FILE__),\
          std::to_string(__LINE__)).ostream()

#define OPENMVG_LOG_IF(MODE, condition) \
  ((MODE) < openMVG::system::logger::logger_severity || !(condition))? (void) 0\
  : openMVG::system::logger::NullBuffer() &\
      openMVG::system::logger::StreamMessageLogger().setMessage(\
        MODE,\
        openMVG::system::logger::filename(__FILE__),\
        std::to_string(__LINE__)).ostream()


#define OPENMVG_LOG_INFO    OPENMVG_LOG(openMVG::system::logger::ELogMode::VERBOSITY_INFO)
#define OPENMVG_LOG_WARNING OPENMVG_LOG(openMVG::system::logger::ELogMode::VERBOSITY_WARNING)
#define OPENMVG_LOG_ERROR   OPENMVG_LOG(openMVG::system::logger::ELogMode::VERBOSITY_ERROR)

#define OPENMVG_LOG_INFO_IF(condition)    OPENMVG_LOG_IF(openMVG::system::logger::ELogMode::VERBOSITY_INFO, condition)
#define OPENMVG_LOG_WARNING_IF(condition) OPENMVG_LOG_IF(openMVG::system::logger::ELogMode::VERBOSITY_WARNING, condition)
#define OPENMVG_LOG_ERROR_IF(condition)   OPENMVG_LOG_IF(openMVG::system::logger::ELogMode::VERBOSITY_ERROR, condition)

#endif // OPENMVG_SYSTEM_LOGGER_HPP
