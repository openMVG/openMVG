
// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "software/VO/AlternativeVO/VOFolderProcessor.hpp"
#include "software/VO/Abstract_Tracker.hpp"
#include "software/VO/Tracker.hpp"
#if defined HAVE_OPENCV
#include "software/VO/Tracker_opencv_klt.hpp"
#endif
#include "openMVG/image/image_io.hpp"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

namespace alternative_vo
{
VOFolderProcessor::VOFolderProcessor( const std::string & inputFolder )
  : m_input_folder( inputFolder ) ,
    m_current_file_ID( 0 )
{
  m_input_files = stlplus::folder_files( m_input_folder );
  // clean invalid image file
  {
    std::vector<std::string> vec_image_;
    for ( size_t i = 0; i < m_input_files.size(); ++i )
    {
      if ( openMVG::image::GetFormat( m_input_files[i].c_str() ) != openMVG::image::Unknown )
      {
        vec_image_.push_back( m_input_files[i] );
      }
    }
    vec_image_.swap( m_input_files );
  }
  std::sort( m_input_files.begin(), m_input_files.end() );

  // Initialize the Visual Odometry interface
  m_tracker.reset( new openMVG::VO::Tracker_fast_dipole );
  m_monocular_vo.reset( new openMVG::VO::VO_Monocular( m_tracker.get(), 1500 ) );
}

/**
* @brief Reset processing to the begining
*/
void VOFolderProcessor::Reset( void )
{
  m_current_file_ID = 0;
}

/**
* @brief get current filename beeing processed
*/
std::string VOFolderProcessor::CurrentFileName( void )
{
  return m_input_files[ m_current_file_ID ];
}

/**
* @brief get full filename for a given id
* @retval empty string if id is out of range
* @retval full path for given id if id is in valid range
*/
std::string VOFolderProcessor::FullFileName( const size_t id )
{
  if (id >= m_input_files.size() )
  {
    return "";
  }
  else
  {
    return stlplus::create_filespec( m_input_folder, m_input_files[ id ] );
  }
}


/**
* @brief Get total number of frame
*/
size_t VOFolderProcessor::NbFrame( void )
{
  return m_input_files.size();
}


/**
 * @brief get current frame ID
 */
size_t VOFolderProcessor::CurrentFrameID( void )
{
  return m_current_file_ID;
}


/**
* @brief Try to step forward
* @retval false if process has ended
* @retval true if new image could be processed
*/
bool VOFolderProcessor::StepForward( void )
{
  if (m_current_file_ID + 1 >= m_input_files.size() )
  {
    return false;
  }
  else
  {
    const std::string sImageFilename = stlplus::create_filespec( m_input_folder, m_input_files[ m_current_file_ID ] );

    openMVG::image::Image<unsigned char> currentImage;

    if ( openMVG::image::ReadImage( sImageFilename.c_str(), &currentImage ) )
    {
      m_monocular_vo->nextFrame( currentImage , m_current_file_ID );
    }

    ++m_current_file_ID;
    return true;
  }
}

/**
* @brief Get current position of the tracked points
*/
std::vector< VOViewerPoint > VOFolderProcessor::GetCurrentTrackedPoints( void ) const
{
  std::vector< VOViewerPoint > res;
  res.reserve(m_monocular_vo->landmark_.size());

  for ( size_t idx = 0; idx < m_monocular_vo->landmark_.size(); ++idx )
  {
    if ( std::find( m_monocular_vo->trackedLandmarkIds_.begin(),
                    m_monocular_vo->trackedLandmarkIds_.end(), idx )
         == m_monocular_vo->trackedLandmarkIds_.end() )
    {
      continue;
    }

    const openMVG::VO::Landmark & landmark = m_monocular_vo->landmark_[idx];
    if ( landmark.obs_.back().frameId_ == m_current_file_ID - 1 && landmark.obs_.size() > 1 )
    {
      const std::deque<openMVG::VO::Measurement> & obs = landmark.obs_;

      // draw the current tracked point
      {
        std::deque<openMVG::VO::Measurement>::const_reverse_iterator iter = obs.rbegin();
        VOViewerPoint cur_tracked;
        cur_tracked.m_color = VOViewerPoint::DEFAULT_TRACKED_POINT_COLOR;
        cur_tracked.m_pt = iter->pos_;
        res.push_back( cur_tracked );
      }
    }
  }
  return res;
}

/**
* @brief Get current position of the newly created points in the last processed frame
*/
std::vector< VOViewerPoint > VOFolderProcessor::GetCreatedPoints( void ) const
{
  std::vector< VOViewerPoint > res;
  res.reserve(m_monocular_vo->landmark_.size());

  for ( size_t idx = 0; idx < m_monocular_vo->landmark_.size(); ++idx )
  {
    if ( std::find( m_monocular_vo->trackedLandmarkIds_.begin(),
                    m_monocular_vo->trackedLandmarkIds_.end(), idx )
         == m_monocular_vo->trackedLandmarkIds_.end() )
    {
      continue;
    }

    const openMVG::VO::Landmark & landmark = m_monocular_vo->landmark_[idx];
    if ( landmark.obs_.back().frameId_ == m_current_file_ID - 1 && landmark.obs_.size() > 1 )
    {

    }
    else // Draw the new initialized point
    {
      if ( landmark.obs_.size() == 1 )
      {
        const std::deque<openMVG::VO::Measurement> & obs = landmark.obs_;
        std::deque<openMVG::VO::Measurement>::const_iterator iter = obs.begin();

        VOViewerPoint cur_created;
        cur_created.m_color = VOViewerPoint::DEFAULT_NEW_POINT_COLOR;
        cur_created.m_pt = iter->pos_;
        res.push_back( cur_created );
      }
    }
  }
  return res;
}

/**
* @brief Get for each tracked points, their full trajectories
*/
std::vector< VOViewerLine > VOFolderProcessor::GetCurrentTrackTrajectories( void ) const
{
  std::vector< VOViewerLine > res;
  res.reserve(m_monocular_vo->landmark_.size());

  for ( size_t idx = 0; idx < m_monocular_vo->landmark_.size(); ++idx )
  {
    if ( std::find( m_monocular_vo->trackedLandmarkIds_.begin(),
                    m_monocular_vo->trackedLandmarkIds_.end(), idx )
         == m_monocular_vo->trackedLandmarkIds_.end() )
    {
      continue;
    }

    const openMVG::VO::Landmark & landmark = m_monocular_vo->landmark_[idx];
    if ( landmark.obs_.back().frameId_ == m_current_file_ID - 1 && landmark.obs_.size() > 1 )
    {
      const std::deque<openMVG::VO::Measurement> & obs = landmark.obs_;

      std::deque<openMVG::VO::Measurement>::const_reverse_iterator iter = obs.rbegin();
      std::deque<openMVG::VO::Measurement>::const_reverse_iterator iterEnd = obs.rend();

      VOViewerLine cur_traj;
      cur_traj.m_color = VOViewerLine::DEFAULT_LINE_COLOR;
      /* Trajectories */
      int limit = 10;
      for (; iter != iterEnd && limit >= 0; ++iter, --limit )
      {
        const openMVG::Vec2f & p0 = iter->pos_;
        cur_traj.m_pts.push_back( p0 );
      }
      res.push_back( cur_traj );
    }
  }
  return res;
}


}
