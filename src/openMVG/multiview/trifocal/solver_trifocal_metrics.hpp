// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
//:\file
//\author Ricardo Fabbri Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Apr 19 20:25:06 -03 2022
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Pierre MOULON
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRIFOCAL_METRICS_HPP
#define OPENMVG_MULTIVIEW_TRIFOCAL_METRICS_HPP

#include "openMVG/system/logger.hpp"
#include "openMVG/multiview/trifocal/trifocal_model.hpp"

namespace openMVG {
namespace trifocal {
  
struct NormalizedSquaredPointReprojectionOntoOneViewError {
  static double Error(
    const trifocal_model_t &tt,
    const Vec &bearing_0, // x,y,tangentialx,tangentialy
    const Vec &bearing_1,
    const Vec &bearing_2);

  // Meant to be run by the 3 points given to trifocal solver
  static bool  Check(
    const trifocal_model_t &tt,
    const Vec &bearing_0, // x,y,tangentialx,tangentialy
    const Vec &bearing_1,
    const Vec &bearing_2);

  // get a reasonable error threshold in normalized coordinates
  //
  // take a (threshold,0) vector along the x axis and 
  // transform to normalized coordinates. Currently ignores skew
  // 
  // todo(better guess is possible or use angular error)
  inline static double threshold_pixel_to_normalized(double threshold, const double k[2][3]) {
    return threshold/k[0][0];
  }

  inline static double threshold_normalized_to_pixel(double threshold, const double k[2][3]) {
    return threshold*k[0][0];
  }
};


} // namespace trifocal
} // namespace OpenMVG

#endif  // OPENMVG_MULTIVIEW_TRIFOCAL_METRICS_HPP
