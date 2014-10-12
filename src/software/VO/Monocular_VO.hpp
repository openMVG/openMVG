
// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef MONOCULAR_VO_HPP
#define MONOCULAR_VO_HPP

#include <deque>

namespace openMVG  {
namespace VO  {

/// Store an image observation of a landmark
struct Measurement
{
  Measurement(size_t frameId, const Vec2f & p): _frameId(frameId), _pos(p) { }

  size_t _frameId;
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
template <typename Tracker_T>
struct VO_Monocular
{
  std::deque<Landmark> _landmark;
  std::vector<size_t> _trackedLandmarkIds;

  Tracker_T _tracker;
  int _maxTrackedFeatures ;
  std::vector<PointFeature> _pt_to_track, _pt_tracked;
  std::vector<unsigned char> _tracking_status;

  VO_Monocular(const int maxTrackedFeatures = 500) : _maxTrackedFeatures(maxTrackedFeatures)
  {
  }

  bool nextFrame(const Image<unsigned char> & ima, const size_t frameId)
  {
    bool bTrackerStatus = _tracker.track(ima, _pt_to_track, _pt_tracked, _tracking_status);

    {
      // Update landmark observation
      // Update tracking point set (if necessary)
      // Compute pose if ?
      // Refine if ?


      //-- Update landmark observation
      for(unsigned int i = 0; i < _tracking_status.size(); ++i)
      {
        if (_tracking_status[i])
        {
          // if tracked
          const PointFeature & a = _pt_to_track[i];
          const PointFeature & b = _pt_tracked[i];
          // enable a spatial filter
          if (DistanceL2(a.coords(), b.coords()) < 50)
          {
            const size_t tracker_landmark_id = _trackedLandmarkIds[i];
            _landmark[tracker_landmark_id]._obs.push_back(Measurement(frameId, b.coords()));
            _pt_to_track[i] = b; // update the tracked point
          }
          else // too much motion
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
      const size_t countTracked = accumulate(_tracking_status.begin(), _tracking_status.end(), 0);
      std::cout << "#tracked: " << countTracked << std::endl;

      // Update tracking point set (if necessary)
      if (countTracked < _maxTrackedFeatures)
      {
        if (_pt_to_track.empty()) {
          _pt_to_track.resize(_maxTrackedFeatures);
          _pt_tracked.resize(_maxTrackedFeatures);
          _trackedLandmarkIds.resize(_maxTrackedFeatures);
          _tracking_status.resize(_maxTrackedFeatures);
          std::fill(_tracking_status.begin(), _tracking_status.end(), (unsigned char)0);
        }

        // add some new feature
        const size_t count = _maxTrackedFeatures - countTracked;
        std::vector<PointFeature> new_pt;
        _tracker.detect(ima, new_pt, count);
        std::cout << "#features added: " << new_pt.size() << std::endl;
        size_t j = 0;
        for(unsigned int i = 0; i < _tracking_status.size(); ++i)
        {
          if (!_tracking_status[i])
          {
            // Create a new landmark
            Landmark landmark;
            landmark._obs.push_back(Measurement(frameId, new_pt[j].coords()));
            _landmark.push_back(landmark);
            // a new landmark ID have be tracked
            _trackedLandmarkIds[i] = _landmark.size()-1;

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
  }
};

} // namespace VO
} // namespace openMVG


#endif // MONOCULAR_VO_HPP
