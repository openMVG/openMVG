// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_BA_HPP
#define OPENMVG_SFM_DATA_BA_HPP

#include "openMVG/cameras/Camera_Common.hpp"

namespace openMVG {
namespace sfm {

class SfM_Data;

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

/// Structure to control which parameter will be refined during the BundleAjdustment process
struct Optimize_Options
{
  cameras::Intrinsic_Parameter_Type intrinsics;
  Extrinsic_Parameter_Type extrinsics;
  Structure_Parameter_Type structure;

  Optimize_Options
  (
    cameras::Intrinsic_Parameter_Type intrinsics_opt = cameras::Intrinsic_Parameter_Type::ADJUST_ALL,
    Extrinsic_Parameter_Type extrinsics_opt = Extrinsic_Parameter_Type::ADJUST_ALL,
    Structure_Parameter_Type structure_opt = Structure_Parameter_Type::ADJUST_ALL
  )
  :intrinsics(intrinsics_opt),
   extrinsics(extrinsics_opt),
   structure(structure_opt)
  {
  }
};

class Bundle_Adjustment
{
  public:
  // Perform a Bundle Adjustment on the SfM scene (refinement only asked parameters)
  virtual bool Adjust
  (
    // the SfM scene to refine
    SfM_Data & sfm_data,
    // tell which parameter needs to be adjusted
    const Optimize_Options options
  ) = 0;

  // TODO: Use filter to say wich parameter is const or not (allow to refine only a subpart of the intrinsics or the poses)
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_BA_HPP
