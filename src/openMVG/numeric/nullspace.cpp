// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2017 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/numeric/nullspace.hpp"

namespace openMVG
{

double Nullspace
(
  const Eigen::Ref<const Mat> & A,
  Eigen::Ref<Vec> nullspace
)
{
  if ( A.rows() >= A.cols() )
  {
    Eigen::JacobiSVD<Mat> svd( A, Eigen::ComputeFullV );
    nullspace = svd.matrixV().col( A.cols() - 1 );
    return svd.singularValues()( A.cols() - 1 );
  }
  // Extend A with rows of zeros to make it square. It's a hack, but it is
  // necessary until Eigen supports SVD with more columns than rows.
  Mat A_extended( A.cols(), A.cols() );
  A_extended.block( A.rows(), 0, A.cols() - A.rows(), A.cols() ).setZero();
  A_extended.block( 0, 0, A.rows(), A.cols() ) = A;
  return Nullspace( A_extended, nullspace );
}

double Nullspace2
(
  const Eigen::Ref<const Mat> & A,
  Eigen::Ref<Vec> x1,
  Eigen::Ref<Vec> x2
)
{
  if ( A.rows() >= A.cols() )
  {
    Eigen::JacobiSVD<Mat> svd( A, Eigen::ComputeFullV );
    const Mat & V = svd.matrixV();
    x1 = V.col( A.cols() - 1 );
    x2 = V.col( A.cols() - 2 );
    return svd.singularValues()( A.cols() - 1 );
  }
  // Extend A with rows of zeros to make it square. It's a hack, but it is
  // necessary until Eigen supports SVD with more columns than rows.
  Mat A_extended( A.cols(), A.cols() );
  A_extended.block( A.rows(), 0, A.cols() - A.rows(), A.cols() ).setZero();
  A_extended.block( 0, 0, A.rows(), A.cols() ) = A;
  return Nullspace2( A_extended, x1, x2 );
}


} // namespace openMVG

