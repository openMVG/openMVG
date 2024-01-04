// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MATCHING_KVLD_DRAW_H
#define OPENMVG_MATCHING_KVLD_DRAW_H

#include <algorithm>
#include <vector>

#include "openMVG/image/image_drawing.hpp"
#include "openMVG/features/feature.hpp"

namespace openMVG {

//-- A slow but accurate way to draw K-VLD lines
void getKVLDMask(
  image::Image<unsigned char> *maskL,
  image::Image<unsigned char> *maskR,
  const std::vector<features::SIOPointFeature > &vec_F1,
  const std::vector<features::SIOPointFeature > &vec_F2,
  const std::vector<Pair>& vec_matches,
  const std::vector<bool>& vec_valide,
  const openMVG::Mat& mat_E)
{
  for (size_t it1 = 0; it1 < vec_matches.size() - 1; it1++ )
  {
    for (size_t it2 = it1 + 1; it2 < vec_matches.size(); it2++ )
    {
      if (vec_valide[ it1 ] && vec_valide[ it2 ] && mat_E( it1, it2 ) >= 0 )
      {
        const features::SIOPointFeature & l1 = vec_F1[ vec_matches[ it1 ].first ];
        const features::SIOPointFeature & l2 = vec_F1[ vec_matches[ it2 ].first ];
        const float l = ( l1.coords() - l2.coords() ).norm();
        const int widthL = std::max( 1.f, l / ( dimension + 1.f ) );

        image::DrawLineThickness(l1.x(), l1.y(), l2.x(), l2.y(), 255, widthL, maskL);

        const features::SIOPointFeature & r1 = vec_F2[ vec_matches[ it1 ].second ];
        const features::SIOPointFeature & r2 = vec_F2[ vec_matches[ it2 ].second ];
        const float r = ( r1.coords() - r2.coords() ).norm();
        const int widthR = std::max( 1.f, r / ( dimension + 1.f ) );

        image::DrawLineThickness(r1.x(), r1.y(), r2.x(), r2.y(), 255, widthR, maskR);
      }
    }
  }
}

}; // namespace openMVG

#endif // OPENMVG_MATCHING_KVLD_DRAW_H
