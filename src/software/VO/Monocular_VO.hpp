
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef MONOCULAR_VO_HPP
#define MONOCULAR_VO_HPP

#include "openMVG/features/features.hpp"
#include <openMVG/numeric/numeric.h>

#include <software/VO/Abstract_Tracker.hpp>
#include <deque>
#include <set>

namespace openMVG  {
namespace VO  {

/// Store an image observation of a landmark
struct Measurement
{
  Measurement
  (
    const uint32_t & frameId,
    const Vec2f & p
  ): _frameId(frameId), _pos(p)
  { }
  Measurement( const Measurement & src ) = default ;

  uint32_t _frameId;
  Vec2f _pos;
};

/// A 3D point with it's associated image observations
struct Landmark
{
  Landmark():_pt(-1,-1,-1) {}

  Vec3 _pt;
  std::deque<Measurement> _obs;
};

/// Monocular test interface

struct VO_Monocular
{
  // Structure and visibility
  std::deque<Landmark> _landmark;
  std::vector<uint64_t> _trackedLandmarkIds;

  // Landmark Id per frame (for easier pairing)
  std::deque< std::set<uint32_t> > _landmarkListPerFrame;

  // Tracking
  Abstract_Tracker * tracker_;
  int _maxTrackedFeatures ;
  std::vector<features::PointFeature> _pt_to_track, _pt_tracked;
  std::vector<bool> _tracking_status;

  VO_Monocular
  (
    Abstract_Tracker * tracker,
    const int maxTrackedFeatures = 1500
    // Add an abstract camera model here
  )
  : _maxTrackedFeatures(maxTrackedFeatures),
    tracker_(tracker)
  {
  }

  bool nextFrame
  (
    const image::Image<unsigned char> & ima,
    const size_t frameId
  )
  {
    const bool bTrackerStatus = tracker_->track(ima, _pt_to_track, _pt_tracked, _tracking_status);
    std::cout << (int) bTrackerStatus  << " : tracker status" << std::endl;
    _landmarkListPerFrame.push_back(std::set<uint32_t>());
    if (_landmarkListPerFrame.size()==1 || bTrackerStatus)
    {
      // Update landmark observation
      // Update tracking point set (if necessary)
      // Compute pose if ?
      // Refine if ?

      //-- Update landmark observation
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (int i = 0; i < (int)_tracking_status.size(); ++i)
      {
        if (_tracking_status[i])
        {
          // if tracked
          const features::PointFeature & a = _pt_to_track[i];
          const features::PointFeature b = _pt_tracked[i];
          // Use a spatial filter
          if ((a.coords() - b.coords()).norm() < 25)
          {
            const size_t tracker_landmark_id = _trackedLandmarkIds[i];
            #ifdef OPENMVG_USE_OPENMP
            #pragma omp critical
            #endif
            {
              _landmark[tracker_landmark_id]._obs.push_back(Measurement( frameId , b.coords() ) );
              _landmarkListPerFrame.back().insert(tracker_landmark_id);
            }
            _pt_to_track[i] = b; // update the tracked point
          }
          else // the feature does not longer appear or tracking failed to find the feature
          {
            _tracking_status[i] = 0;
          }
        }
        else
        {
          // feature was not tracked
        }
      }

      // Count the number of tracked features
      const size_t countTracked = std::accumulate(_tracking_status.begin(), _tracking_status.end(), 0);
      std::cout << "#tracked: " << countTracked << std::endl;

      // try compute pose and decide if it's a Keyframe
      if (frameId > 0 && _landmarkListPerFrame.size() > 1)
      {
        size_t lastKf = frameId-1;
        std::vector<size_t> ids;
        std::set_intersection(
          _landmarkListPerFrame[lastKf].begin(), _landmarkListPerFrame[lastKf].end(),
          _landmarkListPerFrame[frameId].begin(), _landmarkListPerFrame[frameId].end(),
          std::back_inserter(ids));
        std::cout << "Track in common with the last Keyframe: " << ids.size() << std::endl;
      }

      // Update tracking point set (if necessary)
      if (countTracked < _maxTrackedFeatures)
      {
        if (_pt_to_track.empty())
        {
          _pt_to_track.resize(_maxTrackedFeatures);
          _pt_tracked.resize(_maxTrackedFeatures);
          _trackedLandmarkIds.resize(_maxTrackedFeatures);
          _tracking_status.resize(_maxTrackedFeatures);
          std::fill(_tracking_status.begin(), _tracking_status.end(), false);
        }

        // add some new feature
        const size_t count = _maxTrackedFeatures - countTracked;
        std::vector<features::PointFeature> new_pt;
        tracker_->detect(ima, new_pt, count);
        std::cout << "#features added: " << new_pt.size() << std::endl;
        size_t j = 0;
        for (size_t i = 0; i < _tracking_status.size(); ++i)
        {
          if (!_tracking_status[i])
          {
            // Create a new landmark
            Landmark landmark;
            landmark._obs.emplace_back(frameId, new_pt[j].coords());
            _landmark.push_back(landmark);
            // a new landmark ID have be tracked

            _trackedLandmarkIds[i] = _landmark.size()-1;
            _landmarkListPerFrame.back().insert(_landmark.size() - 1);

            _pt_to_track[i] = new_pt[j];
            ++j;
          }
        }
        std::cout << "_landmark.size() " << _landmark.size() << std::endl;
      }
    }
    if (!bTrackerStatus)
    {
       // re-localization
    }
    return bTrackerStatus;
  }
};

} // namespace VO
} // namespace openMVG


#endif // MONOCULAR_VO_HPP
