// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 nomoko AG, Sebastien Chappuis<sebastien@nomoko.camera>, Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/types.hpp"
#include "openMVG/sfm/sfm_data_BA_ceres.hpp"// we only need this one for IntrinsicsToCostFunction and for ceres options
#include <set>

namespace ceres {
  class Problem;
}

namespace openMVG{
namespace sfm{

// Normal bundle adjustment but with separator landmarks
// positions fixed !
class Bundle_Adjustment_Fixed_Separators
{
  protected:
    Bundle_Adjustment_Ceres::BA_Ceres_options ceres_options_;

  public:
  Bundle_Adjustment_Fixed_Separators(const Bundle_Adjustment_Ceres::BA_Ceres_options & options = Bundle_Adjustment_Ceres::BA_Ceres_options());

  Bundle_Adjustment_Ceres::BA_Ceres_options & ceres_options(){return ceres_options_;}

  bool Adjust
  (
   SfM_Data & sfm_data, // scene to refine
   const std::set<IndexT> & separator_tracks_ids
  );

protected:
  virtual void configureProblem(
      ceres::Problem & problem,
      SfM_Data & sfm_data,
      const std::set<IndexT>& separator_tracks_ids,
      Hash_Map<IndexT, std::vector<double> > & map_intrinsics,
      Hash_Map<IndexT, std::vector<double> > & map_poses);
};


}
}
