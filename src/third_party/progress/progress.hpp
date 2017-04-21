
/**
 * @file progress.h
 * @brief Interface for progress update
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

#include <atomic>
#include <string>
#include <algorithm>

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
    ~C_Progress() = default;

    /** @brief Initializer of the C_Progress class
     * @param expected_count The number of step of the process
     * @param msg an optional status message
     **/
    virtual void restart ( unsigned long ulExpected_count, const std::string& msg=std::string())
    //  Effects: display appropriate scale
    //  Postconditions: count()==0, expected_count()==expected_count
    {
      _count = _next_tic_count = _tic = 0;
      _expected_count = ulExpected_count;

      if ( !_expected_count ) _expected_count = 1;  // prevent divide by zero
    } // restart

    /** @brief Indicator if the current operation should be aborted.
     * @return Return true if the process has been canceled by the user.
     **/
    virtual bool hasBeenCanceled()const { return false; }
    /**
     * @brief Post-Increment operator
     * @param increment the number of step that we want to increment the internal step counter
     * @return the value of the internal count => count()
     **/
    unsigned long operator+= ( unsigned long ulIncrement )
    //  Effects: Increment appropriate progress tic if needed.
    //  Postconditions: count()== original count() + increment
    //  Returns: count().
    {
      if ( ( _count += ulIncrement ) >= _next_tic_count ) { inc_tic(); }
      return _count;
    }

    /** @brief A dummy progress reporter. Does nothing.
     **/
    static C_Progress& dummy()
    {
        static C_Progress null;
        return null;
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
        long unsigned int iModulusStep = std::max(_expected_count / iModulus, (long unsigned int)1);  //1 prevent divide per 0
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
    std::atomic<unsigned long> _count, _expected_count, _next_tic_count;
    std::atomic<unsigned int> _tic;

  private:
    virtual void inc_tic()
    {
      _next_tic_count = static_cast<unsigned long> ( ( _tic/50.0 ) *_expected_count );
    } // inc_tic
};

#endif
