
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_ROTATION_AVERAGING_H_
#define OPENMVG_MULTIVIEW_ROTATION_AVERAGING_H_

#include "openMVG/numeric/numeric.h"
#include <vector>
#include <map>

//--
//-- Implementation related to rotation averaging.
// . Compute global rotation from a list of relative estimates.
//
//- Implementation of algorithm from Thesis titled:
//- [1] "Robust Multiview Reconstruction."
//- Author : Daniel Martinec.
//- Date : July 2, 2008.
//--
namespace openMVG   {
namespace rotation_averaging  {

// [1] 6.7.2 Consistent Rotation page 89
// Closest Rotation Estimation R = U*transpose(V)
//  approximate rotation in the Frobenius norm using SVD
inline Mat3 ClosestSVDRotationMatrix(const Mat3 & rotMat)
{
  // Closest orthogonal matrix
  Eigen::JacobiSVD<Mat3> svd(rotMat,Eigen::ComputeFullV|Eigen::ComputeFullU);
  Mat3 U = svd.matrixU();
  Mat3 V = svd.matrixV();
  return U*V.transpose();
}

typedef std::vector<std::pair<std::pair<size_t, size_t>, Mat3> > vector_RelativeRotMotion;

//-- Solve the Global Rotation matrix registration for each camera given a list
//    of relative orientation using matrix parametrization
//    [1] formula 6.62 page 100. Dense formulation.
//- nCamera:               The number of camera to solve
//- vec_rotationEstimate:  The relative rotation i->j
//- vec_ApprRotMatrix:     The output global rotation

// Minimization of the norm of:
// => || rj - Rij * ri ||= 0
// With rj et rj the global rotation and Rij the relative rotation from i to j.
//
// Example:
// 0_______2
//  \     /
//   \   /
//    \ /
//     1
//
// nCamera = 3
// vector.add( make_pair(make_pair(0,1), R01) );
// vector.add( make_pair(make_pair(1,2), R12) );
// vector.add( make_pair(make_pair(0,2), R02) );
//
static bool L2RotationAveraging( size_t nCamera,
  const vector_RelativeRotMotion & vec_relativeRotEstimate,
  std::vector<Mat3> & vec_ApprRotMatrix)
{
  using namespace std;
  Mat3 Id = Mat3::Identity();
  const size_t nRotationEstimation = vec_relativeRotEstimate.size();
  //--
  // Setup the Action Matrix
  //--
  // nCamera * 3 because each columns have 3 elements.
  Mat A = Mat::Zero(nRotationEstimation*3, 3*nCamera);
  //-- Encode constraint (6.62 Martinec Thesis page 100):
  size_t cpt = 0;
  for(vector_RelativeRotMotion::const_iterator
    iter = vec_relativeRotEstimate.begin();
    iter != vec_relativeRotEstimate.end();
    iter++, cpt++)
  {
    const std::pair<std::pair<size_t, size_t>, Mat3> & Elem = *iter;

    //-- Encode rj - Rij * ri = 0
    int i = Elem.first.first;
    int j = Elem.first.second;

    A.block<3,3>(3 * cpt, 3 * i) = - Elem.second;
    A.block<3,3>(3 * cpt, 3 * j) =   Id;
  }

  // Solve Ax=0 => SVD
  Eigen::JacobiSVD<Mat> svd(A,Eigen::ComputeFullV);
  const Vec & NullspaceVector0 = svd.matrixV().col(A.cols()-1);
  const Vec & NullspaceVector1 = svd.matrixV().col(A.cols()-2);
  const Vec & NullspaceVector2 = svd.matrixV().col(A.cols()-3);

  //--
  // Search the closest matrix :
  //  - From solution of SVD get back column and reconstruct Rotation matrix
  //  - Enforce the orthogonality constraint
  //     (approximate rotation in the Frobenius norm using SVD).
  //--
  vec_ApprRotMatrix.clear();
  vec_ApprRotMatrix.reserve(nCamera);
  for(size_t i=0; i < nCamera; ++i)
  {
    Mat3 Rotation;
    Rotation << NullspaceVector0.segment(3 * i, 3),
                NullspaceVector1.segment(3 * i, 3),
                NullspaceVector2.segment(3 * i, 3);

    //-- Compute the closest SVD rotation matrix
    Rotation = ClosestSVDRotationMatrix(Rotation);
    vec_ApprRotMatrix.push_back(Rotation);
  }
  // Force R0 to be Identity
  Mat3 R0T = vec_ApprRotMatrix[0].transpose();
  for(size_t i = 0; i < nCamera; ++i) {
    vec_ApprRotMatrix[i] *= R0T;
  }

  return true;
}

} // namespace rotation_averaging
} // namespace openMVG

#endif //OPENMVG_MULTIVIEW_ROTATION_AVERAGING_H_
