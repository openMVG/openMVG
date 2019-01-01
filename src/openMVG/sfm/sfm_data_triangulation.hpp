// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_TRIANGULATION_HPP
#define OPENMVG_SFM_SFM_DATA_TRIANGULATION_HPP

#include <set>

#include "openMVG/numeric/eigen_alias_definition.hpp"
#include "openMVG/sfm/sfm_landmark.hpp"
#include "openMVG/types.hpp"

namespace openMVG { namespace sfm { struct SfM_Data; } }

namespace openMVG {
namespace sfm {

/// Generic basis struct for triangulation of track data contained
///  in the SfM_Data scene structure.
struct SfM_Data_Structure_Computation_Basis
{
  bool bConsole_verbose_;

  explicit SfM_Data_Structure_Computation_Basis(bool bConsoleVerbose = false);

  virtual void triangulate(SfM_Data & sfm_data) const = 0;
};


/// Triangulation of track data contained in the structure of a SfM_Data scene.
// Use a blind estimation:
// - Triangulate tracks using all observations
// - Inlier/Outlier classification is done by a cheirality test
struct SfM_Data_Structure_Computation_Blind: public SfM_Data_Structure_Computation_Basis
{
  explicit SfM_Data_Structure_Computation_Blind(bool bConsoleVerbose = false);

  void triangulate(SfM_Data & sfm_data) const override;
};

/// Triangulation of track data contained in the structure of a SfM_Data scene.
// Use a robust estimation:
// - Triangulate tracks using a RANSAC scheme
// - Check cheirality and a pixel residual error
struct SfM_Data_Structure_Computation_Robust: public SfM_Data_Structure_Computation_Basis
{
  explicit SfM_Data_Structure_Computation_Robust
  (
    const double max_reprojection_error = 4.0, // pixels
    const IndexT min_required_inliers = 3,
    const IndexT min_sample_index = 3,
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
    Landmark & landmark) const;

private:

  // -- DATA
  double max_reprojection_error_;
  const IndexT min_required_inliers_;
  const IndexT min_sample_index_;

};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_TRIANGULATION_HPP
