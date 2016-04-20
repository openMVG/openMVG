// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_TRIANGULATION_HPP
#define OPENMVG_SFM_DATA_TRIANGULATION_HPP

#include "openMVG/sfm/sfm_data.hpp"

namespace openMVG {
namespace sfm {

/// Generic basis struct for triangulation of track data contained
///  in the SfM_Data scene structure.
struct SfM_Data_Structure_Computation_Basis
{
  bool bConsole_verbose_;

  SfM_Data_Structure_Computation_Basis(bool bConsoleVerbose = false);

  virtual void triangulate(SfM_Data & sfm_data) const = 0;
};


/// Triangulation of track data contained in the structure of a SfM_Data scene.
// Use a blind estimation:
// - Triangulate tracks using all observations
// - Inlier/Outlier classification is done by a cheirality test
struct SfM_Data_Structure_Computation_Blind: public SfM_Data_Structure_Computation_Basis
{
  SfM_Data_Structure_Computation_Blind(bool bConsoleVerbose = false);

  void triangulate(SfM_Data & sfm_data) const override;
};

/// Triangulation of track data contained in the structure of a SfM_Data scene.
// Use a robust estimation:
// - Triangulate tracks using a RANSAC scheme
// - Check cheirality and a pixel residual error (TODO: make it a parameter)
struct SfM_Data_Structure_Computation_Robust: public SfM_Data_Structure_Computation_Basis
{
  SfM_Data_Structure_Computation_Robust
  (
    const double max_reprojection_error = 4, // pixels
    bool bConsoleVerbose = false
  );

  void triangulate(SfM_Data & sfm_data) const override;

  /// Robust triangulation of track data contained in the structure
  /// All observations must have View with valid Intrinsic and Pose data
  /// Invalid landmark are removed.
  void robust_triangulation(SfM_Data & sfm_data) const;

  /// Robustly try to estimate the best 3D point using a ransac Scheme
  /// Return true for a successful triangulation
  bool robust_triangulation(
    const SfM_Data & sfm_data,
    const Observations & obs,
    Landmark & landmark,
    const IndexT min_required_inliers = 3,
    const IndexT min_sample_index = 3) const;

private:
  /// Triangulate a given track from a selection of observations
  Vec3 track_sample_triangulation(
    const SfM_Data & sfm_data,
    const Observations & obs,
    const std::set<IndexT> & samples) const;

  // -- DATA
  double max_reprojection_error_;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_TRIANGULATION_HPP
