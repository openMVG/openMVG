// Copyright (c) 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/image_features/image_features.hpp"
#include "openMVG/image_features/detector_fast.hpp"
#include "openMVG/image_features/descriptor_dipole.hpp"

namespace openMVG{
namespace image_features{

class Fast_dipole_img_features{

  int _size, _threshold;

public:

  Fast_dipole_img_features(int size = 9, int threshold = 30):
    _size(size), _threshold(threshold)
  {
  }

  void detect_and_describe(
    const Image<unsigned char> & ima,
    Regions *& regions)
  {
    //a. detect
    //b. describe
    //c. save data to regions container

    //a. detect
    std::vector<SIOPointFeature> feats;
    FastCornerDetector fastCornerDetector(_size, _threshold);
    fastCornerDetector.detect(ima, feats);

    if (regions) delete regions;
    regions = new DISSOCIATED_DIPOLES_Regions;
    DISSOCIATED_DIPOLES_Regions * regionsCasted = dynamic_cast<DISSOCIATED_DIPOLES_Regions*>(regions);
    regionsCasted->features().resize(feats.size());
    regionsCasted->descriptors().resize(feats.size());

#ifdef USE_OPENMP
  #pragma omp parallel for
#endif
    for (int i = 0; i < static_cast<int>(feats.size()); ++i)
    {
      const SIOPointFeature & feat = feats[i];
      //b. describe (compute the descriptor)
      Descriptor<float, 20> desc;
      PickNaiveDipole(ima, feat.x(), feat.y(), 3.0f, 0.0f, &desc[0]);
      //c. save data to regions container
      regionsCasted->features()[i] = feat;
      regionsCasted->descriptors()[i] = desc;
    }
  }
};

} // image_features
} // namespace openMVG
