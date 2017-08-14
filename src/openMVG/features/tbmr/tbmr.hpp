// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2014-2016 Yongchao Xu, Pascal Monasse, Thierry Géraud, Laurent Najman
// Copyright (c) 2016 Pierre Moulon.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Implements [1]: the *Tree-based Morse Regions* (TBMR) local
// feature detector.
// This detector extracts as features the connected components of the level
// sets of the input intensity image. Among all such regions, the ones that
// are maximally topological equivalent to the critical regions (based on
// Morse theory) are selected.
// TBMRs are affine co-variant, as well as largely co-variant to generic
// diffeomorphic transformations, and truly invariant to affine contrast
// changes.

//------------------
//-- Bibliography --
//------------------
//- [1] "Tree-Based Morse Regions: A Topological Approach to Local Feature Detection"
//- Authors: Yongchao Xu, Pascal Monasse, Thierry Géraud, Laurent Najman
//- Date: 2014
//- IEEE Transactions on Image Processing, Institute of Electrical and Electronics Engineers (IEEE)

#ifndef OPENMVG_FEATURES_TBMR_TBMR_HPP
#define OPENMVG_FEATURES_TBMR_TBMR_HPP

#include <functional>
#include <vector>
#include "openMVG/features/feature.hpp"

namespace openMVG { namespace image { template <typename T> class Image; } }

namespace openMVG
{
namespace features
{
namespace tbmr
{
  /**
  * @brief Extract the TBMR (Tree-Based Morse Regions) of an image
  * @param[in] ima Preceding direction
  * @param[out] features detected TBMR as affine features
  * @param[in] cmp considered pixel ordering
  * @param[in] minimumSize minimum area of a region to be detected (i.e. 30 pixels)
  * @param[in] maximumRelativeSize minimum area of a region to be detected relative to image area
  * @param cmp ordering (std::less => BRIGHT features; or std::greater => DARK features)
  */
  template <typename Ordering = std::less<unsigned char>>
  void Extract_tbmr
  (
    const image::Image<unsigned char> & ima,
    std::vector<features::AffinePointFeature> & features,
    Ordering cmp = Ordering(),
    const unsigned int minimumSize = 30,
    const double maximumRelativeSize = 0.01
  );

} // namespace tbmr
} // namespace features
} // namespace openMVG

#endif // OPENMVG_FEATURES_TBMR_TBMR_HPP
