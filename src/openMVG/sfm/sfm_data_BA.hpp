// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_BA_HPP
#define OPENMVG_SFM_DATA_BA_HPP

namespace openMVG {
namespace sfm {

class Bundle_Adjustment
{
  public:
  // Perform a Bundle Adjustment on the SfM scene with refinement of the requested parameters
  virtual bool Adjust(
    SfM_Data & sfm_data,            // the SfM scene to refine
    bool bRefineRotations = true,   // tell if pose rotations will be refined
    bool bRefineTranslations = true,// tell if the pose translation will be refined
    bool bRefineIntrinsics = true,  // tell if the camera intrinsic will be refined
    bool bRefineStructure = true)   // tell if the structure will be refined
  = 0;

  // TODO: Use filter to say wich parameter is const or not (allow to refine only a subpart of the intrinsics or the poses)
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_BA_HPP
