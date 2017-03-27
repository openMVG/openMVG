// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (C) 2013  "Robot Vision Group, NLPR, CASIA", Zhenhua Wang,
// Bin Fan and Fuchao Wu.

// Copyright (C) 2014: Adaptation to openMVG by Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_FEATURES_LIOP_LIOP_DESCRIPTOR_HPP
#define OPENMVG_FEATURES_LIOP_LIOP_DESCRIPTOR_HPP


//------------------
//-- Bibliography --
//------------------
//- [1] "Local Intensity Order Pattern for Feature Description"
//- Authors: Zhenhua Wang, Bin Fan and Fuchao Wu
//- Date: 2011, ICCV, IEEE International Conference on Computer Vision


namespace openMVG { namespace features { class SIOPointFeature; } }
namespace openMVG { namespace image { template <typename T> class Image; } }

#include <map>
#include <vector>

namespace openMVG {
namespace features{
namespace LIOP    {

class Liop_Descriptor_Extractor
{
private:
  std::map<int, unsigned char> m_LiopPatternMap;
  std::vector<int> m_LiopPosWeight;
public:

  Liop_Descriptor_Extractor();

  void extract(
    const image::Image<unsigned char> & I,
    const SIOPointFeature & feat,
    float desc[144]);

  void CreateLIOP_GOrder(
    const image::Image<float> & outPatch,
    const image::Image<unsigned char> & flagPatch,
    const int inRadius,
    float desc[144]) const;

  void GeneratePatternMap(
    std::map<int,unsigned char> & pattern_map,
    std::vector<int> & pos_weight,
    unsigned char n);
};

} // namespace LIOP
} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_LIOP_LIOP_DESCRIPTOR_HPP
