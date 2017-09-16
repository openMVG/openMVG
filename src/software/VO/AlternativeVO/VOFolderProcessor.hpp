
// Copyright (c) 2016 Romuald Perrot

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ALTERNATIVE_VO_VO_FOLDER_PROCESSOR_HPP_
#define ALTERNATIVE_VO_VO_FOLDER_PROCESSOR_HPP_

#include "software/VO/AlternativeVO/VOViewerDrawableElements.hpp"

#include "software/VO/Monocular_VO.hpp"

#include <string>
#include <vector>

namespace openMVG  {
namespace VO  {

struct Abstract_Tracker;

}
}

namespace alternative_vo
{
/**
* @brief Class encapsulating VO processing of an image folder
*/
class VOFolderProcessor
{
  public:

    VOFolderProcessor( const std::string & inputFolder );

    /**
    * @brief Reset processing to the begining
    */
    void Reset( void );

    /**
    * @brief get current filename beeing processed
    */
    std::string CurrentFileName( void );

    /**
    * @brief get full filename for a given id
    * @retval empty string if id is out of range
    * @retval full path for given id if id is in valid range
    */
    std::string FullFileName( const size_t id );

    /**
    * @brief get current frame ID
    */
    size_t CurrentFrameID( void );

    /**
    * @brief Get total number of frame
    */
    size_t NbFrame( void );

    /**
    * @brief Try to step forward
    * @retval false if process has ended
    * @retval true if new image could be processed
    */
    bool StepForward( void );

    /**
    * @brief Get current position of the tracked points
    */
    std::vector< VOViewerPoint > GetCurrentTrackedPoints( void ) const;

    /**
    * @brief Get current position of the newly created points in the last processed frame
    */
    std::vector< VOViewerPoint > GetCreatedPoints( void ) const;

    /**
    * @brief Get for each tracked points, their full trajectories
    */
    std::vector< VOViewerLine > GetCurrentTrackTrajectories( void ) const;

  private:

    std::string m_input_folder;
    std::vector< std::string > m_input_files;

    std::unique_ptr<openMVG::VO::VO_Monocular> m_monocular_vo;
    std::unique_ptr<openMVG::VO::Abstract_Tracker> m_tracker;

    size_t m_current_file_ID;
};
}

#endif
