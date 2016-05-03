// Copyright (c) 2015 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_DATA_BA_HPP
#define OPENMVG_SFM_DATA_BA_HPP

#include <openMVG/stl/bitmask.hpp>

namespace openMVG {
namespace sfm {

enum BA_Refine
{
  BA_REFINE_NONE = 0,
  BA_REFINE_ROTATION = 1, //< refine pose rotations
  BA_REFINE_TRANSLATION = 2, //< refine pose translations
  BA_REFINE_INTRINSICS = 4, //< refine camera intrinsics
  BA_REFINE_STRUCTURE = 8, //<  refine structure (i.e. 3D points)
  /// Refine all parameters
  BA_REFINE_ALL = BA_REFINE_ROTATION | BA_REFINE_TRANSLATION | BA_REFINE_INTRINSICS | BA_REFINE_STRUCTURE,
};

OPENMVG_BITMASK(BA_Refine)

class Bundle_Adjustment
{
  public:
  /**
   * @brief Perform a Bundle Adjustment on the SfM scene with refinement of the requested parameters
   * 
   * @param[in,out] sfm_data: sfm data scene to modify with the BA
   * @param[in] refineOptions: choose what you want to refine
   * @return 
   */
  virtual bool Adjust(
    SfM_Data & sfm_data,
    BA_Refine refineOptions = BA_REFINE_ALL) = 0;

  // TODO: Use filter to say wich parameter is const or not (allow to refine only a subpart of the intrinsics or the poses)
};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_SFM_DATA_BA_HPP
