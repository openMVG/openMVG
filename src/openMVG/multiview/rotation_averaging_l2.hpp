
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_ROTATION_AVERAGING_L2_H_
#define OPENMVG_MULTIVIEW_ROTATION_AVERAGING_L2_H_

#include "openMVG/multiview/rotation_averaging_common.hpp"
#include <vector>
#include <map>

#ifdef _MSC_VER
#pragma warning( once : 4267 ) //warning C4267: 'argument' : conversion from 'size_t' to 'const int', possible loss of data
#endif

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
namespace l2  {

// <eigenvalue, eigenvector> pair comparator
bool compare_first_abs(std::pair<double, Vec> const &x, std::pair<double, Vec> const &y)
{
 return fabs(x.first) < fabs(y.first);
}

// [1] 6.7.2 Consistent Rotation page 89
// Closest Rotation Estimation R = U*transpose(V)
//  approximate rotation in the Frobenius norm using SVD
inline Mat3 ClosestSVDRotationMatrix(const Mat3 & rotMat)
{
  // Closest orthogonal matrix
  Eigen::JacobiSVD<Mat3> svd(rotMat,Eigen::ComputeFullV|Eigen::ComputeFullU);
  const Mat3 U = svd.matrixU();
  const Mat3 V = svd.matrixV();
  return U*V.transpose();
}

//-- Solve the Global Rotation matrix registration for each camera given a list
//    of relative orientation using matrix parametrization
//    [1] formula 6.62 page 100. Dense formulation.
//- nCamera:               The number of camera to solve
//- vec_rotationEstimate:  The relative rotation i->j
//- vec_ApprRotMatrix:     The output global rotation

// Minimization of the norm of:
// => || wij * (rj - Rij * ri) ||= 0
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
// vector.add( RelRotationData(0,1, R01) );
// vector.add( RelRotationData(1,2, R12) );
// vector.add( RelRotationData(0,2, R02) );
//
static bool L2RotationAveraging( size_t nCamera,
  const std::vector<RelRotationData>& vec_relativeRot,
  // Output
  std::vector<Mat3> & vec_ApprRotMatrix)
{
  const size_t nRotationEstimation = vec_relativeRot.size();
  //--
  // Setup the Action Matrix
  //--
  std::vector<Eigen::Triplet<double> > tripletList;
  tripletList.reserve(nRotationEstimation*12); // 3*3 + 3
  //-- Encode constraint (6.62 Martinec Thesis page 100):
  sMat::Index cpt = 0;
  for(std::vector<RelRotationData>::const_iterator
    iter = vec_relativeRot.begin();
    iter != vec_relativeRot.end();
    iter++, cpt++)
  {
   const RelRotationData & Elem = *iter;

   //-- Encode weight * ( rj - Rij * ri ) = 0
   const size_t i = iter->i;
   const size_t j = iter->j;

   // A.block<3,3>(3 * cpt, 3 * i) = - Rij * weight;
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt, 3 * i, - iter->Rij(0,0) * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt, 3 * i + 1, - iter->Rij(0,1) * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt, 3 * i + 2, - iter->Rij(0,2) * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt + 1, 3 * i, - iter->Rij(1,0) * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt + 1, 3 * i + 1, - iter->Rij(1,1) * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt + 1, 3 * i + 2, - iter->Rij(1,2) * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt + 2, 3 * i, - iter->Rij(2,0) * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt + 2, 3 * i + 1, - iter->Rij(2,1) * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt + 2, 3 * i + 2, - iter->Rij(2,2) * iter->weight));

   // A.block<3,3>(3 * cpt, 3 * j) = Id * weight;
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt, 3 * j, 1.0 * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt + 1, 3 * j + 1, 1.0 * iter->weight));
   tripletList.push_back(Eigen::Triplet<double>(3 * cpt + 2, 3 * j + 2, 1.0 * iter->weight));
  }

  // nCamera * 3 because each columns have 3 elements.
  sMat A(nRotationEstimation*3, 3*nCamera);
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  tripletList.clear();

  sMat AtAsparse = A.transpose() * A;
  const Mat AtA = Mat(AtAsparse); // convert to dense

  // You can use either SVD or eigen solver (eigen solver will be faster) to solve Ax=0

  // Solve Ax=0 => SVD
  //Eigen::JacobiSVD<Mat> svd(A,Eigen::ComputeFullV);
  //const Vec & NullspaceVector0 = svd.matrixV().col(A.cols()-1);
  //const Vec & NullspaceVector1 = svd.matrixV().col(A.cols()-2);
  //const Vec & NullspaceVector2 = svd.matrixV().col(A.cols()-3);

  // Solve Ax=0 => eigen vectors
  Eigen::SelfAdjointEigenSolver<Mat> es(AtA, Eigen::ComputeEigenvectors);

  if (es.info() != Eigen::Success)
  {
    return false;
  }
  else
  {
    // Sort abs(eigenvalues)
    std::vector<std::pair<double, Vec> > eigs(AtA.cols());
    for (size_t i = 0; i < AtA.cols(); ++i)
    {
      eigs[i] = std::make_pair(es.eigenvalues()[i], es.eigenvectors().col(i));
    }
    std::stable_sort(eigs.begin(), eigs.end(), &compare_first_abs);

    const Vec & NullspaceVector0 = eigs[0].second;
    const Vec & NullspaceVector1 = eigs[1].second;
    const Vec & NullspaceVector2 = eigs[2].second;

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
    const Mat3 R0T = vec_ApprRotMatrix[0].transpose();
    for(size_t i = 0; i < nCamera; ++i) {
      vec_ApprRotMatrix[i] *= R0T;
    }

    return true;
  }
}

} // namespace l2
} // namespace rotation_averaging
} // namespace openMVG

#endif //OPENMVG_MULTIVIEW_ROTATION_AVERAGING_L2_H_

