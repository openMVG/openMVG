
/**
 * @file progress_display.h
 * @brief Simple progress bar for console application
 * @author Pierre MOULON
 *
 * Copyright (c) 2011, 2012, 2013 Pierre MOULON
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef PROGRESS_DISPLAY
#define PROGRESS_DISPLAY

#include "third_party/progress/progress.hpp"

#include <iostream>
#include <mutex>

 /**
 * @brief C_Progress_display displays an appropriate indication of progress at an appropriate place in an appropriate form.
 * -- Class modeled from BOOST::PROGRESS::progress_display
 * Usage :
 * C_Progress_display my_progress_bar( fileList.size() );
 * for (list<string>::const_iterator it = fileList.begin(); it != fileList.end(); ++it, ++my_progress_bar)
 * the C_Progress_display::operator++ take in charge the display of the progression (***)
 * {
 *   const string & filename = *it;
 *   ... // do something
 * }
 * The displayed result will be like the following :
 * 0%   10   20   30   40   50   60   70   80   90   100%
 * |----|----|----|----|----|----|----|----|----|----|
 * ************************** ...
 *
 **/

class C_Progress_display : public C_Progress
{
public:
    C_Progress_display(): m_os(std::cout), m_msg("\n")
    {
    }
    /** @brief Standard constructor
    * @param expected_count The number of step of the process
    * @param os the stream that will be used for display
    * @param msg the status string
    **/
    explicit C_Progress_display(
      unsigned long ulExpected_count,
      std::ostream & os = std::cout,
      const std::string & msg = "\n")
      // os is hint; implementation may ignore, particularly in embedded systems
      : m_os(os), m_msg(msg)
    {
      restart(ulExpected_count, msg);
    }

    /** @brief Initializer of the C_Progress_display class
    * @param expected_count The number of step of the process
    * @param msg updates the status string. Can be empty to keep the last one.
    **/
    void restart(unsigned long ulExpected_count, const std::string& msg=std::string()) override
    //  Effects: display appropriate scale
    //  Postconditions: count()==0, expected_count()==expected_count
    {
      C_Progress::restart(ulExpected_count, msg); //-- Initialize the base class
      if (!msg.empty())
        m_msg = msg;
      m_os
        << m_msg
        << "0%   10   20   30   40   50   60   70   80   90   100%\n"
        <<  "|----|----|----|----|----|----|----|----|----|----|"
        << std::endl;  // endl implies flush, which ensures display
    } // restart

private:
    /// Internal reference over the stream used to display the progress bar
    std::ostream &     m_os;  // may not be present in all imps
    // String used to formatting display
    std::string  m_msg;
    // Mutex to protect concurrent writting in the console
    std::mutex m_mutex_display;

    /** @brief Function that check if we have to append an * in the progress bar function
    **/
    void inc_tic() override
    {
      std::lock_guard<std::mutex> console_display_lock(m_mutex_display);
      // use of floating point ensures that both large and small counts
      // work correctly.
      unsigned int tics_needed = static_cast<unsigned int> ((static_cast<double> (_count) /_expected_count) *50.0);
      do { m_os << '*' << std::flush; } while (++_tic < tics_needed);
      _next_tic_count = static_cast<unsigned long> ((_tic/50.0) *_expected_count);
      if (_count == _expected_count)
      {
        if (_tic < 51)
          m_os << '*';
        m_os << std::endl;
      }
    } // display_tic
};

/** @} */

#endif
