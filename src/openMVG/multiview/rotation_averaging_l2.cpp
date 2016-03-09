
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/multiview/rotation_averaging_l2.hpp"
#include <vector>
#include <map>

#include <ceres/ceres.h>
#include <ceres/rotation.h>

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

// [1] 6.7.2 Consistent Rotation page 89
// Closest Rotation Estimation R = U*transpose(V)
//  approximate rotation in the Frobenius norm using SVD
Mat3 ClosestSVDRotationMatrix
(
  const Mat3 & rotMat
)
{
  // Closest orthogonal matrix
  Eigen::JacobiSVD<Mat3> svd(rotMat,Eigen::ComputeFullV|Eigen::ComputeFullU);
  const Mat3 & U = svd.matrixU();
  const Mat3 & V = svd.matrixV();

  Mat3 rot_mat = U*V.transpose();
  // Ensure that the projection will have a +1 determinant (SO3).
  if (rot_mat.determinant() < 0) {
    rot_mat *= -1.0;
  }
  return rot_mat;
}

// <eigenvalue, eigenvector> pair comparator
bool compare_first_abs(std::pair<double, Vec> const &x, std::pair<double, Vec> const &y)
{
 return fabs(x.first) < fabs(y.first);
}

//-- Solve the Global Rotation matrix registration for each camera given a list
//    of relative orientation using matrix parametrization
//    [1] formula 6.62 page 100. Dense formulation.
//- nCamera:               The number of camera to solve
//- vec_rotationEstimate:  The relative rotation i->j
//- vec_ApprRotMatrix:     The output global rotation
//
// Example:
// 0_______2
//  \     /
//   \   /
//    \ /
//     1
//
// nCamera = 3
// vector.add( RelativeRotation(0,1, R01) );
// vector.add( RelativeRotation(1,2, R12) );
// vector.add( RelativeRotation(0,2, R02) );
//
// It solves the following linear system:
// => wij * (R_j - R_{i,j} * R_i) = 0
// Note:
// R_j, R_j are the global rotations
// R_{i,j} are the relative rotations from i to j.
//
// It creates a (3 * m) x (3 * n) system.
// => m = #R_{i,j} => the number of relative rotations.
// => n => the number of view (camera)
//
bool L2RotationAveraging
(
  size_t nCamera,
  const RelativeRotations& vec_relativeRot,
  // Output
  std::vector<Mat3> & global_rotations
)
{
  const size_t nRotationEstimation = vec_relativeRot.size();
  //--
  // Setup the Action Matrix
  //--
  std::vector<Eigen::Triplet<double> > tripletList;
  tripletList.reserve(nRotationEstimation*12); // 3*3 + 3
  //-- Encode constraint (6.62 Martinec Thesis page 100):
  sMat::Index cpt = 0;
  for(RelativeRotations::const_iterator
    iter = vec_relativeRot.begin();
    iter != vec_relativeRot.end();
    iter++, cpt++)
  {
    const RelativeRotation & Elem = *iter;

    //-- Encode weight * ( rj - Rij * ri ) = 0
    const sMat::Index i = iter->i;
    const sMat::Index j = iter->j;

    // A.block<3,3>(3 * cpt, 3 * i) = - Rij * weight;
    tripletList.emplace_back(3 * cpt, 3 * i, - iter->Rij(0,0) * iter->weight);
    tripletList.emplace_back(3 * cpt, 3 * i + 1, - iter->Rij(0,1) * iter->weight);
    tripletList.emplace_back(3 * cpt, 3 * i + 2, - iter->Rij(0,2) * iter->weight);
    tripletList.emplace_back(3 * cpt + 1, 3 * i, - iter->Rij(1,0) * iter->weight);
    tripletList.emplace_back(3 * cpt + 1, 3 * i + 1, - iter->Rij(1,1) * iter->weight);
    tripletList.emplace_back(3 * cpt + 1, 3 * i + 2, - iter->Rij(1,2) * iter->weight);
    tripletList.emplace_back(3 * cpt + 2, 3 * i, - iter->Rij(2,0) * iter->weight);
    tripletList.emplace_back(3 * cpt + 2, 3 * i + 1, - iter->Rij(2,1) * iter->weight);
    tripletList.emplace_back(3 * cpt + 2, 3 * i + 2, - iter->Rij(2,2) * iter->weight);

   // A.block<3,3>(3 * cpt, 3 * j) = Id * weight;
    tripletList.emplace_back(3 * cpt, 3 * j, iter->weight);
    tripletList.emplace_back(3 * cpt + 1, 3 * j + 1, iter->weight);
    tripletList.emplace_back(3 * cpt + 2, 3 * j + 2, iter->weight);
  }

  // nCamera * 3 because each columns have 3 elements.
  Mat AtA(3*nCamera,3*nCamera);
  {
    sMat A(nRotationEstimation*3, 3*nCamera);
    A.setFromTriplets(tripletList.begin(), tripletList.end());
    tripletList.clear();
    tripletList.shrink_to_fit();

    const sMat AtAsparse = A.transpose() * A;
    AtA = Mat(AtAsparse); // convert to dense
  }

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
    global_rotations.clear();
    global_rotations.reserve(nCamera);
    for(size_t i=0; i < nCamera; ++i)
    {
      Mat3 Rotation;
      Rotation << NullspaceVector0.segment(3 * i, 3),
                  NullspaceVector1.segment(3 * i, 3),
                  NullspaceVector2.segment(3 * i, 3);

      //-- Compute the closest SVD rotation matrix
      global_rotations.emplace_back(ClosestSVDRotationMatrix(Rotation));
    }
    // Force R0 to be Identity
    const Mat3 R0T = global_rotations[0].transpose();
    for(size_t i = 0; i < nCamera; ++i) {
      global_rotations[i] *= R0T;
    }
  }
  return true;
}

// Ceres Functor to minimize global rotation regarding fixed relative rotation
struct CeresPairRotationError {
  CeresPairRotationError(const openMVG::Vec3& relative_rotation,  const double weight)
    :relative_rotation_(relative_rotation), weight_(weight) {}

  // The error is given by the rotation cycle error (R2 * R1.t) * RRel.t
  template <typename T>
  bool operator()
  (
    const T* angleAxis1,
    const T* angleAxis2,
    T* residuals
  )
  const
  {
    const T relative_rotation[3] = {
      T(relative_rotation_[0]),
      T(relative_rotation_[1]),
      T(relative_rotation_[2])};

    // Convert angle axis rotations to rotation matrices
    Eigen::Matrix<T, 3, 3> RRel, R1, R2;
    ceres::AngleAxisToRotationMatrix( relative_rotation, RRel.data());
    ceres::AngleAxisToRotationMatrix( angleAxis1, R1.data());
    ceres::AngleAxisToRotationMatrix( angleAxis2, R2.data());

    // Compute the "cycle" rotation error for the given edge:
    //  relative error between two given global rotations and the relative one
    const Eigen::Matrix<T, 3, 3> cycle_rotation_mat = (R2 * R1.transpose()) * RRel.transpose();
    Eigen::Matrix<T, 3, 1> cycle_rotation;
    ceres::RotationMatrixToAngleAxis( cycle_rotation_mat.data(), cycle_rotation.data());

    residuals[0] = T(weight_) * cycle_rotation(0);
    residuals[1] = T(weight_) * cycle_rotation(1);
    residuals[2] = T(weight_) * cycle_rotation(2);

    return true;
  }

  const openMVG::Vec3 relative_rotation_;
  const double weight_;
};

bool L2RotationAveraging_Refine
(
  const RelativeRotations & vec_relativeRot,
  std::vector<openMVG::Mat3> & vec_ApprRotMatrix
)
{
  if (vec_relativeRot.size() == 0 ||vec_ApprRotMatrix.size() == 0 ) {
    std::cout << "Skip nonlinear rotation optimization, no sufficient data provided " << std::endl;
    return false;
  }

  // Convert global rotation to AngleAxis representation
  std::vector<openMVG::Vec3> vec_Rot_AngleAxis(vec_ApprRotMatrix.size());
  for (int i=0; i < vec_ApprRotMatrix.size(); ++i)
  {
    ceres::RotationMatrixToAngleAxis((const double*)vec_ApprRotMatrix[i].data(), vec_Rot_AngleAxis[i].data());
  }

  // Setup the problem and a robust loss function.
  ceres::Problem problem;
  const double robust_loss_width = 0.03; // 2° error along one axis (perhaps a bit too strict)

  for (size_t ii = 0; ii < vec_relativeRot.size(); ++ii)
  {
    const int i = vec_relativeRot[ii].i;
    const int j = vec_relativeRot[ii].j;
    const openMVG::Mat3& Rrel = vec_relativeRot[ii].Rij;
    const double w = vec_relativeRot[ii].weight;

    openMVG::Vec3 rotAngleAxis;
    ceres::RotationMatrixToAngleAxis((const double*)Rrel.data(), rotAngleAxis.data());

    ceres::CostFunction* cost_function =
      new ceres::AutoDiffCostFunction<CeresPairRotationError, 3, 3, 3>(
      new CeresPairRotationError(rotAngleAxis, w));

    ceres::LossFunction* loss_function = new ceres::SoftLOneLoss(w * robust_loss_width);
    problem.AddResidualBlock(cost_function,
      loss_function,
      vec_Rot_AngleAxis[i].data(),
      vec_Rot_AngleAxis[j].data());
  }
  ceres::Solver::Options solverOptions;
  solverOptions.minimizer_progress_to_stdout = false;
  solverOptions.logging_type = ceres::SILENT;
  // Since the problem is sparse, use a sparse solver iff available
  if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
  {
    solverOptions.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
    solverOptions.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  }
  else if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
  {
    solverOptions.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
    solverOptions.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  }
  else if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::EIGEN_SPARSE))
  {
    solverOptions.sparse_linear_algebra_library_type = ceres::EIGEN_SPARSE;
    solverOptions.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
  }
  else
  {
    solverOptions.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
  }
#ifdef OPENMVG_USE_OPENMP
  solverOptions.num_threads = omp_get_max_threads();
  solverOptions.num_linear_solver_threads = omp_get_max_threads();
#endif // OPENMVG_USE_OPENMP

  ceres::Solver::Summary summary;
  ceres::Solve(solverOptions, &problem, &summary);
  // std::cout << summary.FullReport() << std::endl;

  if (summary.IsSolutionUsable())
  {
    // Convert back the AngleAxis rotations to rotations matrices
    for (int i=0; i < vec_ApprRotMatrix.size(); ++i)
    {
      ceres::AngleAxisToRotationMatrix(
        vec_Rot_AngleAxis[i].data(), vec_ApprRotMatrix[i].data());
    }
  }
  return summary.IsSolutionUsable();
}

} // namespace l2
} // namespace rotation_averaging
} // namespace openMVG


