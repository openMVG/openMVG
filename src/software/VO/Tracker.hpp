// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef TRACKER_VO_HPP
#define TRACKER_VO_HPP
#include <software/VO/Abstract_Tracker.hpp>
#include <openMVG/features/dipole/dipole_descriptor.hpp>

#include <openMVG/features/fast/fast_detector.hpp>
#include <openMVG/features/feature.hpp>
#include <openMVG/features/feature_container.hpp>
#include "openMVG/matching/metric.hpp"

#include <random>
#include <vector>

namespace openMVG  {
namespace VO  {

// Implement tracking by description:
//  - Each tracked point uses a local description and is tracked thanks to descriptor matching
struct Tracker_fast_dipole : public Abstract_Tracker
{
  // data for tracking
  image::Image<unsigned char> _prev_img;

  /// Try to track current point set in the provided image
  /// return false when tracking failed (=> to send frame to relocalization)
  bool track
  (
    const image::Image<unsigned char> & ima,
    const std::vector<features::PointFeature> & pt_to_track,
    std::vector<features::PointFeature> & pt_tracked,
    std::vector<bool> & status
  ) override
  {
    if (!pt_to_track.empty())
    {
      const features::PointFeatures & _prevPts = pt_to_track;

      //-- Compute descriptors for the previous tracked point and perform matching
      std::vector<float> prev_descriptors(20*pt_to_track.size());
      for (size_t i=0; i < pt_to_track.size(); ++i)
      {
        features::PickASDipole(_prev_img, _prevPts[i].x(), _prevPts[i].y(), 10.5f, 0.0f, &prev_descriptors[i*20]);
      }

      features::PointFeatures current_feats;
      features::FastCornerDetector fastCornerDetector(9, 5);
      fastCornerDetector.detect(ima, current_feats);
      std::vector<float> current_descriptors(20*current_feats.size());
      for (size_t i=0; i < current_feats.size(); ++i)
      {
        features::PickASDipole(ima, current_feats[i].x(), current_feats[i].y(), 10.5f, 0.0f, &current_descriptors[i*20]);
      }

      // Compute the matches
      {
        #ifdef OPENMVG_USE_OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (int i=0; i < (int)pt_to_track.size(); ++i)
        {
          size_t best_idx = std::numeric_limits<size_t>::infinity();

          typedef openMVG::matching::L2<float> metricT;
          metricT metric;
          metricT::ResultType best_distance = 30;//std::numeric_limits<double>::infinity();
          for (size_t j=0; j < current_feats.size(); ++j)
          {
            // Spatial filter
            //if ((pt_to_track[i].coords() - current_feats[j].coords()).norm() > 50)
            //  continue; // Too much distance between the feat, it is not necessary to compute the descriptor distance

            metricT::ResultType distance = metric(&prev_descriptors[i*20], &current_descriptors[j*20], 20);
            if (distance < best_distance)
            {
              best_idx = j;
              best_distance = distance;
            }
          }
          if (best_idx != std::numeric_limits<size_t>::infinity())
          {
            pt_tracked[i].coords() << current_feats[best_idx].x(), current_feats[best_idx].y();
            status[i] = true;
          }
          else
          {
            status[i] = false;
          }
        }
      }
    }
    // swap frame for the next tracking iteration
    _prev_img = ima;

    const size_t tracked_point_count = std::accumulate(status.begin(), status.end(), 0);
    return (tracked_point_count != 0);
  }

  // suggest new feature point for tracking (count point are kept)
  bool detect
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::PointFeature> & pt_to_track,
    const size_t count
  ) const override
  {
    std::mt19937 gen(std::mt19937::default_seed);
    const std::vector<unsigned char> scores = {30, 20, 10, 5};
    // use a sequence of scores 'in order to deal with lighting change'
    features::PointFeatures feats;
    for (const unsigned char fast_score : scores)
    {
      features::FastCornerDetector fastCornerDetector(9, fast_score);
      fastCornerDetector.detect(ima, feats);

      if (feats.size() > count)
      {
        // shuffle to avoid to sample only in one bucket
        std::shuffle(feats.begin(), feats.end(), gen);
        feats.resize(count); // cut the array to keep only a given count of features
        pt_to_track.swap(feats);

        return true;
      }
      feats.clear();
    }
    return false; // Cannot compute a sufficient number of points for the given image
  }
};

} // namespace VO
} // namespace openMVG


#endif // TRACKER_VO_HPP
