// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef MONOCULAR_VO_HPP
#define MONOCULAR_VO_HPP

#include <deque>
#include <set>
#include <numeric>

#include "openMVG/features/feature.hpp"
#include "openMVG/image/image_container.hpp"
#include "openMVG/numeric/eigen_alias_definition.hpp"

#include "software/VO/Abstract_Tracker.hpp"

namespace openMVG  {
namespace VO  {

/// Store an image observation of a landmark
struct Measurement
{
  Measurement
  (
    const uint32_t & frameId,
    const Vec2f & p
  ): frameId_(frameId), pos_(p)
  { }
  Measurement( const Measurement & src ) = default;

  uint32_t frameId_;
  Vec2f pos_;
};

/// A 3D point with it's associated image observations
struct Landmark
{
  Landmark():pt_(-1,-1,-1) {}

  Vec3 pt_;
  std::deque<Measurement> obs_;
};

/// Monocular test interface

struct VO_Monocular
{
  // Structure and visibility
  std::deque<Landmark> landmark_;
  std::vector<uint64_t> trackedLandmarkIds_;

  // Landmark Id per frame (for easier pairing)
  std::deque< std::set<uint32_t> > landmarkListPerFrame_;

  // Tracking
  Abstract_Tracker * tracker_;
  uint32_t maxTrackedFeatures_;
  std::vector<features::PointFeature> pt_to_track_, pt_tracked_;
  std::vector<bool> tracking_status_;

  VO_Monocular
  (
    Abstract_Tracker * tracker,
    const uint32_t maxTrackedFeatures = 1500
    // Add an abstract camera model here
  )
  : tracker_(tracker),
  maxTrackedFeatures_(maxTrackedFeatures)
  {
  }

  bool nextFrame
  (
    const image::Image<unsigned char> & ima,
    const size_t frameId
  )
  {
    const bool bTrackerStatus = tracker_->track(ima, pt_to_track_, pt_tracked_, tracking_status_);
    std::cout << (int) bTrackerStatus  << " : tracker status" << std::endl;
    landmarkListPerFrame_.push_back(std::set<uint32_t>());
    if (landmarkListPerFrame_.size()==1 || bTrackerStatus)
    {
      // Update landmark observation
      // Update tracking point set (if necessary)
      // Compute pose if ?
      // Refine if ?

      //-- Update landmark observation
      #ifdef OPENMVG_USE_OPENMP
      #pragma omp parallel for schedule(dynamic)
      #endif
      for (int i = 0; i < (int)tracking_status_.size(); ++i)
      {
        if (tracking_status_[i])
        {
          // if tracked
          const features::PointFeature & a = pt_to_track_[i];
          const features::PointFeature b = pt_tracked_[i];
          // Use a spatial filter
          if ((a.coords() - b.coords()).norm() < 25)
          {
            const uint64_t tracker_landmark_id = trackedLandmarkIds_[i];
            #ifdef OPENMVG_USE_OPENMP
            #pragma omp critical
            #endif
            {
              landmark_[tracker_landmark_id].obs_.push_back(Measurement( frameId , b.coords() ) );
              landmarkListPerFrame_.back().insert(tracker_landmark_id);
            }
            pt_to_track_[i] = b; // update the tracked point
          }
          else // the feature does not longer appear or tracking failed to find the feature
          {
            tracking_status_[i] = 0;
          }
        }
        else
        {
          // feature was not tracked
        }
      }

      // Count the number of tracked features
      const size_t countTracked = std::accumulate(tracking_status_.begin(), tracking_status_.end(), 0);
      std::cout << "#tracked: " << countTracked << std::endl;

      // try compute pose and decide if it's a Keyframe
      if (frameId > 0 && landmarkListPerFrame_.size() > 1)
      {
        size_t lastKf = frameId-1;
        std::vector<uint32_t> ids;
        std::set_intersection(
          landmarkListPerFrame_[lastKf].begin(), landmarkListPerFrame_[lastKf].end(),
          landmarkListPerFrame_[frameId].begin(), landmarkListPerFrame_[frameId].end(),
          std::back_inserter(ids));
        std::cout << "Track in common with the last Keyframe: " << ids.size() << std::endl;
      }

      // Update tracking point set (if necessary)
      if (countTracked < maxTrackedFeatures_)
      {
        if (pt_to_track_.empty())
        {
          pt_to_track_.resize(maxTrackedFeatures_);
          pt_tracked_.resize(maxTrackedFeatures_);
          trackedLandmarkIds_.resize(maxTrackedFeatures_);
          tracking_status_.resize(maxTrackedFeatures_);
          std::fill(tracking_status_.begin(), tracking_status_.end(), false);
        }

        // add some new feature
        const size_t count = maxTrackedFeatures_ - countTracked;
        std::vector<features::PointFeature> new_pt;
        if (tracker_->detect(ima, new_pt, count))
        {
          std::cout << "#features added: " << new_pt.size() << std::endl;
          size_t j = 0;
          for (size_t i = 0; i < tracking_status_.size(); ++i)
          {
            if (!tracking_status_[i])
            {
              // Create a new landmark
              Landmark landmark;
              landmark.obs_.emplace_back(frameId, new_pt[j].coords());
              landmark_.push_back(landmark);
              // a new landmark ID have be tracked

              trackedLandmarkIds_[i] = landmark_.size()-1;
              landmarkListPerFrame_.back().insert(landmark_.size() - 1);

              pt_to_track_[i] = new_pt[j];
              ++j;
            }
          }
          std::cout << "_landmark.size() " << landmark_.size() << std::endl;
        }
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
