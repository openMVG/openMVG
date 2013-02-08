
/**
 * @file progress.h
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

#ifndef PROGRESS
#define PROGRESS

#include <iostream>
#include <string>
#include <algorithm>
using namespace std;

/**
  * @brief C_Progress manage the appropriate count_step of progress
  *  -- Class modeled from BOOST::PROGRESS::progress_display
  * Usage :
  * C_Progress progress( COUNT );
  * for(int i=0; i < COUNT; ++i, ++progress)
  * {
  *  ... //do something
  * }
  * -- you can access to internal data :
  * //cout<< "Count : " << progress.count() << " pourcent " << progress.count()/(float)progress.expected_count()*100;
 **/
class C_Progress
{
  public:
    /** @brief Standard constructor
     * @param expected_count The number of step of the process
     **/
    explicit C_Progress ( unsigned long ulExpected_count=1 )
    { restart ( ulExpected_count ); }

    ///Destructor
    ~C_Progress() {};

    /** @brief Initializer of the C_Progress class
     * @param expected_count The number of step of the process
     **/
    virtual void           restart ( unsigned long ulExpected_count )
    //  Effects: display appropriate scale
    //  Postconditions: count()==0, expected_count()==expected_count
    {
      _count = _next_tic_count = _tic = 0;
      _expected_count = ulExpected_count;

      if ( !_expected_count ) _expected_count = 1;  // prevent divide by zero
    } // restart

    /**
     * @brief Post-Increment operator
     * @param increment the number of step that we want to increment the internal step counter
     * @return the value of the internal count => count()
     **/
    unsigned long  operator+= ( unsigned long ulIncrement )
    //  Effects: Increment appropriate progress tic if needed.
    //  Postconditions: count()== original count() + increment
    //  Returns: count().
    {
      if ( ( _count += ulIncrement ) >= _next_tic_count ) { inc_tic(); }
      return _count;
    }

    /**
     * @brief Allow to know if the actual pourcentage is the modulus. Used to know when we must display the pourcentage of a progress counter.
     * @param modulus The modulus we want to check.
     * @return Return true if the internal count step is a modulus of the given value or if it is equal to 0.
     **/
    bool isModulus ( int iModulus=10 ) const
    {
      bool bResult=false;
      if(_count==0)
      {
        bResult=true;
      }
      else
      {
        long unsigned int iModulusStep = max(_expected_count / iModulus, (long unsigned int)1);  //1 prevent divide per 0
        bResult=((_count % iModulusStep) == (iModulusStep-1));
      }


      return bResult;
    }

    /**
     * @brief Pre-Increment operator
     * @return the value of _count
     **/
    unsigned long  operator++()           { return operator+= ( 1 ); }

    //-- Accessor
    /**
     * @brief Get the internal count step.
     * @return The pourcentage of the progress => a value in [0;_expected_count].
     **/
    unsigned long  count() const          { return _count; }

    /**
     * @brief Get the _expected_count
     * @return The value of _expected_count
     **/
    unsigned long  expected_count() const { return _expected_count; }

    /**
     * @brief Get the pourcent progress of the process
     * @return The pourcentage of the progress => a value in [0;100].
     **/
    unsigned int  pourcent() const { return (int)( _count/ ( float ) _expected_count*100+.5 ); }

  protected:
    /// Internal data to _count the number of step (the _count can go to the _expected_count value).
    unsigned long _count, _expected_count, _next_tic_count;
    unsigned int  _tic;

  private:
    virtual void inc_tic()
    {
      _next_tic_count = static_cast<unsigned long> ( ( _tic/50.0 ) *_expected_count );
    } // inc_tic
};


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
    /** @brief Standard constructor
      * @param expected_count The number of step of the process
      * @param os the stream that will be used for display
      * @param s1 the starting string
      * @param s2 String used before the |--| sequence
      * @param s3 String used after the  |--| sequence
      **/
    explicit C_Progress_display ( unsigned long ulExpected_count,
                                     std::ostream & os = std::cout,
                                     const std::string & s1 = "\n", //leading strings
                                     const std::string & s2 = "",
                                     const std::string & s3 = "" )
    // os is hint; implementation may ignore, particularly in embedded systems
    : m_os ( os ), m_s1 ( s1 ), m_s2 ( s2 ), m_s3 ( s3 )  { restart ( ulExpected_count ); }

    /** @brief Initializer of the C_Progress_display class
     * @param expected_count The number of step of the process
     **/
    void restart ( unsigned long ulExpected_count )
    //  Effects: display appropriate scale
    //  Postconditions: count()==0, expected_count()==expected_count
    {
      C_Progress::restart ( ulExpected_count ); //-- Initialize the base class

      m_os << m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
      << m_s2 << "|----|----|----|----|----|----|----|----|----|----|"
      << std::endl  // endl implies flush, which ensures display
      << m_s3;
    } // restart

  private:
    /// Internal reference over the stream used to display the progress bar
    std::ostream &     m_os;  // may not be present in all imps
    /// String used to formatting display // string is more general, safer than const char *, and efficiency or size are not issues
    const std::string  m_s1, m_s2, m_s3;

    /** @brief Function that check if we have to append an * in the progress bar function
    **/
    void inc_tic()
    {
      // use of floating point ensures that both large and small counts
      // work correctly.  static_cast<>() is also used several places
      // to suppress spurious compiler warnings.
      unsigned int tics_needed = static_cast<unsigned int> ( ( static_cast<double> ( _count ) /_expected_count ) *50.0 );
      do { m_os << '*' << std::flush; }
      while ( ++_tic < tics_needed );
      _next_tic_count = static_cast<unsigned long> ( ( _tic/50.0 ) *_expected_count );
      if ( _count == _expected_count )
      {
        if ( _tic < 51 )
          m_os << '*';
        m_os << std::endl;
      }
    } // display_tic
};

/** @} */

#endif
