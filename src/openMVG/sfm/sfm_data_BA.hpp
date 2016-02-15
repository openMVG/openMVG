// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_BA_HPP
#define OPENMVG_SFM_DATA_BA_HPP

namespace openMVG {
namespace sfm {

// Enum to control which parameter will be refined during the BA
enum Parameter_Adjustment_Option
{
  ADJUST_CAMERA_ROTATION     = 0x01,
  ADJUST_CAMERA_TRANSLATION  = 0x02,
  ADJUST_CAMERA_INTRINSIC    = 0x04,
  ADJUST_STRUCTURE    = 0x08,
  ADJUST_ALL = ADJUST_CAMERA_ROTATION | ADJUST_CAMERA_TRANSLATION | ADJUST_CAMERA_INTRINSIC | ADJUST_STRUCTURE,
  ADJUST_MOTION_AND_STRUCTURE = ADJUST_CAMERA_ROTATION | ADJUST_CAMERA_TRANSLATION | ADJUST_STRUCTURE
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
    Parameter_Adjustment_Option adjustment_option = ADJUST_ALL
  ) = 0;

  // TODO: Use filter to say wich parameter is const or not (allow to refine only a subpart of the intrinsics or the poses)
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_BA_HPP
