// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2016 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRIANGULATION_NVIEW_HPP
#define OPENMVG_MULTIVIEW_TRIANGULATION_NVIEW_HPP

#include <utility>
#include <vector>

#include "openMVG/numeric/eigen_alias_definition.hpp"

namespace openMVG {

  /// Compute a 3D position of a point from several images of it. In particular,
  ///  compute the projective point X in R^4 such that x = PX.
  /// Algorithm is the standard DLT; for derivation see appendix of Keir's thesis.
  void TriangulateNView
  (
    const Mat3X &x, // x's are landmark bearing vectors in each camera
    const std::vector< Mat34 > &Ps, // Ps are projective cameras
    Vec4 *X
  );

  // This method uses the algebraic distance approximation.
  // Note that this method works better when the 2D points are normalized
  // with an isotopic normalization.
  void TriangulateNViewAlgebraic
  (
    const Mat3X &x, // x's are landmark bearing vectors in each camera
    const std::vector< Mat34 > &Ps, // Ps are projective cameras.
    Vec4 *X
  );

  //Iterated linear method
  class Triangulation
  {
    public:

    size_t size() const;

    void clear();

    void add
    (
      const Mat34& projMatrix,
      const Vec2 & p
    );

    // Return squared L2 sum of error
    double error(const Vec3 &X) const;

    // Compute the corresponding 3D point
    Vec3 compute(int iter = 3) const;

    ////////////////////////////////////////////////////////////////////////////
    // Accessors

    // These values are defined after a successful call to compute
    double minDepth() const { return zmin; }
    double maxDepth() const { return zmax; }
    double error()    const { return err; }

    ////////////////////////////////////////////////////////////////////////////
    // Data members

    protected:
      mutable double zmin; // min depth, mutable since modified in compute(...) const;
      mutable double zmax; // max depth, mutable since modified in compute(...) const;
      mutable double err;  // re-projection error, mutable since modified in compute(...) const;
      std::vector< std::pair<Mat34, Vec2>, Eigen::aligned_allocator<std::pair<Mat34, Vec2>> > views; // Proj matrix and associated image point
  };

}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TRIANGULATION_NVIEW_HPP
