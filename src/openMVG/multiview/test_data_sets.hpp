
// Copyright (c) 2010 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TEST_DATA_SETS_HP
#define OPENMVG_MULTIVIEW_TEST_DATA_SETS_HP

#include <string>
#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {

// A N-view metric dataset.
// All points are seen by all cameras.
struct NViewDataSet {
  std::vector<Mat3> _K;   // Internal parameters (fx, fy, etc).
  std::vector<Mat3> _R;   // Rotation.
  std::vector<Vec3> _t;   // Translation.
  std::vector<Vec3> _C;   // Camera centers.
  Mat3X _X;          // 3D points.
  std::vector<Mat2X> _x;  // Projected points; may have noise added.
  std::vector<Vecu>  _x_ids;// Indexes of points corresponding to the projections

  size_t _n;  // Actual number of cameras.

  //-- Return P=K*[R|t] for the Inth camera
  Mat34 P(size_t i) const;

  /// Export in PLY the point structure and camera and camera looking dir.
  void ExportToPLY(const std::string & out_file_name) const;
};

struct nViewDatasetConfigurator
{
  /// Internal camera parameters (focal, principal point)
  int _fx, _fy, _cx, _cy;

  /// Camera random position parameters
  double _dist;
  double _jitter_amount;

  nViewDatasetConfigurator(int fx = 1000,  int fy = 1000,
                           int cx = 500,   int cy  = 500,
                           double distance = 1.5,
                           double jitter_amount = 0.01 );
};

/// Place cameras on a circle with point in the center
NViewDataSet NRealisticCamerasRing(size_t nviews, size_t npoints,
                                   const nViewDatasetConfigurator
                                     config = nViewDatasetConfigurator());

/// Place cameras on cardiod shape with point in the center
NViewDataSet NRealisticCamerasCardioid(size_t nviews, size_t npoints,
                                       const nViewDatasetConfigurator
                                        config = nViewDatasetConfigurator());

} // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TEST_DATA_SETS_HP
