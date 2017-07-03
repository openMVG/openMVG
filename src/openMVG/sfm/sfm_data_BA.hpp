// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_SFM_DATA_BA_HPP
#define OPENMVG_SFM_SFM_DATA_BA_HPP

#include "openMVG/cameras/Camera_Common.hpp"

namespace openMVG {
namespace sfm {

struct SfM_Data;

/// Enum to control which parameter(s) of the Camera motion must be refined or not
enum class Extrinsic_Parameter_Type : int
{
  NONE                = 0x01,     // Extrinsic parameters will be considered as FIXED
  ADJUST_ROTATION     = 0x02,
  ADJUST_TRANSLATION  = 0x04,
  ADJUST_ALL = ADJUST_ROTATION | ADJUST_TRANSLATION
};

/// Enum to control if the Structure must be refined or not
enum class Structure_Parameter_Type : bool
{
  NONE = false, // Structure will be held as constant
  ADJUST_ALL = true
};

/// Structure to tell to BA if GCP must be use and with which weight
struct Control_Point_Parameter
{
  Control_Point_Parameter
  (
    double weight_val = 20.0,
    bool use_control_points = false
  ): weight(weight_val), bUse_control_points(use_control_points)
  {}
  double weight;
  bool bUse_control_points;
};

/// Structure to control which parameter will be refined during the BundleAjdustment process
struct Optimize_Options
{
  cameras::Intrinsic_Parameter_Type intrinsics_opt;
  Extrinsic_Parameter_Type extrinsics_opt;
  Structure_Parameter_Type structure_opt;
  Control_Point_Parameter control_point_opt;
  bool use_motion_priors_opt;

  Optimize_Options
  (
    const cameras::Intrinsic_Parameter_Type intrinsics = cameras::Intrinsic_Parameter_Type::ADJUST_ALL,
    const Extrinsic_Parameter_Type extrinsics = Extrinsic_Parameter_Type::ADJUST_ALL,
    const Structure_Parameter_Type structure = Structure_Parameter_Type::ADJUST_ALL,
    const Control_Point_Parameter control_point = Control_Point_Parameter(0.0, false), // Default setting does not use GCP in the BA
    const bool use_motion_priors = false
  )
  :intrinsics_opt(intrinsics),
   extrinsics_opt(extrinsics),
   structure_opt(structure),
   control_point_opt(control_point),
   use_motion_priors_opt(use_motion_priors)
  {
  }
};

class Bundle_Adjustment
{
  public:

  virtual ~Bundle_Adjustment() = default;

  // Perform a Bundle Adjustment on the SfM scene (refinement only asked parameters)
  virtual bool Adjust
  (
    // the SfM scene to refine
    sfm::SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options & options
  ) = 0;
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_SFM_DATA_BA_HPP
