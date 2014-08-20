
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_SFM_ROBUST_H
#define OPENMVG_SFM_ROBUST_H

#include "openMVG/numeric/numeric.h"

#include "openMVG/features/feature.hpp"
#include "openMVG/matching/indMatch.hpp"
#include "openMVG/multiview/solver_resection_kernel.hpp"
#include "openMVG/multiview/solver_resection_p3p.hpp"
#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/triangulation.hpp"
#include "openMVG/cameras/PinholeCamera.hpp"

#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

#include "openMVG/bundle_adjustment/problem_data_container.hpp"
#include "openMVG/bundle_adjustment/pinhole_ceres_functor.hpp"

namespace openMVG{
namespace SfMRobust{

using namespace openMVG::matching;
using namespace openMVG::robust;

static const size_t ACRANSAC_ITER = 4096;


/**
 * @brief Estimate the essential matrix from point matches and K matrices.
 *
 * @param[in] K1 camera 1 intrinsics
 * @param[in] K2 camera 2 intrinsics
 * @param[in] x1 camera 1 image points
 * @param[in] x2 camera 2 image points
 * @param[out] pE essential matrix (can be NULL)
 * @param[out] pvec_inliers inliers indices (can be empty)
 * @param[in] size_ima1 width, height of image 1
 * @param[in] size_ima2 width, height of image 2
 * @param[out] errorMax upper bound of the reprojection error of the found solution
 * @param[in] precision upper bound of the desired solution
 */
bool robustEssential(
  const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  Mat3 * pE, std::vector<size_t> * pvec_inliers,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  double * errorMax,
  double precision = std::numeric_limits<double>::infinity())
{
  assert(pvec_inliers != NULL);
  assert(pE != NULL);

  // Use the 5 point solver to estimate E
  typedef openMVG::essential::kernel::FivePointKernel SolverType;
  // Define the AContrario adaptor
  typedef ACKernelAdaptorEssential<
      SolverType,
      openMVG::fundamental::kernel::EpipolarDistanceError,
      UnnormalizerT,
      Mat3>
      KernelType;

  KernelType kernel(x1, size_ima1.first, size_ima1.second,
                    x2, size_ima2.first, size_ima2.second, K1, K2);

  // Robustly estimation of the Essential matrix and it's precision
  std::pair<double,double> acRansacOut = ACRANSAC(kernel, *pvec_inliers,
    ACRANSAC_ITER, pE, precision, false);
  *errorMax = acRansacOut.first;

  return pvec_inliers->size() > 2.5 * SolverType::MINIMUM_SAMPLES;
}

/**
 * @brief Estimate the best possible Rotation/Translation from E.
 *  Four are possible, keep the one with most of the point in front.
 *
 * @param[in] K1 camera 1 intrinsics
 * @param[in] K2 camera 2 intrinsics
 * @param[in] x1 camera 1 image points
 * @param[in] x2 camera 2 image points
 * @param[in] E essential matrix
 * @param[in] vec_inliers inliers indices
 * @param[out] R estimated rotation
 * @param[out] t estimated translation
 */
bool estimate_Rt_fromE(const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  const Mat3 & E, const std::vector<size_t> & vec_inliers,
  Mat3 * R, Vec3 * t)
{
  // Accumulator to find the best solution
  std::vector<size_t> f(4, 0);

  std::vector<Mat3> Es; // Essential,
  std::vector<Mat3> Rs;  // Rotation matrix.
  std::vector<Vec3> ts;  // Translation matrix.

  Es.push_back(E);
  // Recover best rotation and translation from E.
  MotionFromEssential(E, &Rs, &ts);

  //-> Test the 4 solutions will all the point
  assert(Rs.size() == 4);
  assert(ts.size() == 4);

  Mat34 P1, P2;
  Mat3 R1 = Mat3::Identity();
  Vec3 t1 = Vec3::Zero();
  P_From_KRt(K1, R1, t1, &P1);

  for (unsigned int i = 0; i < 4; ++i)
  {
    const Mat3 &R2 = Rs[i];
    const Vec3 &t2 = ts[i];
    P_From_KRt(K2, R2, t2, &P2);
    Vec3 X;

    for (size_t k = 0; k < vec_inliers.size(); ++k)
    {
      const Vec2 & x1_ = x1.col(vec_inliers[k]),
        & x2_ = x2.col(vec_inliers[k]);
      TriangulateDLT(P1, x1_, P2, x2_, &X);
      // Test if point is front to the two cameras.
      if (Depth(R1, t1, X) > 0 && Depth(R2, t2, X) > 0)
      {
          ++f[i];
      }
    }
  }
  // Check the solution:
  const std::vector<size_t>::iterator iter = max_element(f.begin(), f.end());
  if(*iter == 0)
  {
    std::cerr << std::endl << "/!\\There is no right solution,"
      <<" probably intermediate results are not correct or no points"
      <<" in front of both cameras" << std::endl;
    return false;
  }
  const size_t index = std::distance(f.begin(),iter);
  (*R) = Rs[index];
  (*t) = ts[index];

  return true;
}

/// Triangulate a set of points between two view
void triangulate2View_Vector(const Mat34 & P1,
  const Mat34 & P2,
  const std::vector<SIOPointFeature> & vec_feat1,
  const std::vector<SIOPointFeature> & vec_feat2,
  const std::vector<IndMatch>  & vec_index,
  std::vector<Vec3> * pvec_3dPoint,
  std::vector<double> * pvec_residual)
{
  assert(pvec_3dPoint);
  assert(pvec_residual);

  pvec_3dPoint->reserve(vec_index.size());
  pvec_residual->reserve(vec_index.size());

  for (size_t i=0; i < vec_index.size(); ++i)
  {
    //Get corresponding point and triangulate it
    const SIOPointFeature & imaA = vec_feat1[vec_index[i]._i];
    const SIOPointFeature & imaB = vec_feat2[vec_index[i]._j];

    const Vec2 x1 = imaA.coords().cast<double>(),
      x2 = imaB.coords().cast<double>();

    Vec3 X_euclidean = Vec3::Zero();
    TriangulateDLT(P1, x1, P2, x2, &X_euclidean);

    double dResidual2D =
      (PinholeCamera::Residual(P1, X_euclidean, x1) +
      PinholeCamera::Residual(P2, X_euclidean, x2)) /2.0;

    // store 3DPoint and associated residual
    pvec_3dPoint->push_back(X_euclidean);
    pvec_residual->push_back(dResidual2D);
  }

  if (!vec_index.empty())
  {
    double dMin = *min_element(pvec_residual->begin(), pvec_residual->end()),
      dMax = *max_element(pvec_residual->begin(), pvec_residual->end());
    std::cout << std::endl
      << "SfMRobust::triangulate2View_Vector" << std::endl
      << "\t-- Residual min max -- " << dMin <<"\t" << dMax << std::endl;
  }
}

struct ResectionSquaredResidualError {
// Compute the residual of the projection distance(pt2D, Project(P,pt3D))
// Return the squared error
static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D){
  Vec2 x = Project(P, pt3D);
  return (x - pt2D).squaredNorm();
}
};

/// Compute the robust resection of the 3D<->2D correspondences.
bool robustResection(
  const std::pair<size_t,size_t> & imageSize,
  const Mat & pt2D,
  const Mat & pt3D,
  std::vector<size_t> * pvec_inliers,
  const Mat3 * K = NULL,
  Mat34 * P = NULL,
  double * maxError = NULL,
  bool bRefineFocal = true )
{
  double dPrecision = std::numeric_limits<double>::infinity();
  size_t MINIMUM_SAMPLES = 0;
  // Classic resection
  if (K == NULL)
  {
    typedef openMVG::resection::kernel::SixPointResectionSolver SolverType;
    MINIMUM_SAMPLES = SolverType::MINIMUM_SAMPLES;

    typedef ACKernelAdaptorResection<
      SolverType, ResectionSquaredResidualError, UnnormalizerResection, Mat34>
      KernelType;

    KernelType kernel(pt2D, imageSize.first, imageSize.second, pt3D);
    // Robustly estimation of the Projection matrix and it's precision
    std::pair<double,double> ACRansacOut = ACRANSAC(kernel, *pvec_inliers,
      ACRANSAC_ITER, P, dPrecision, true);
    *maxError = ACRansacOut.first;

  }
  else
  {
    //typedef openMVG::euclidean_resection::kernel::EpnpSolver SolverType; // If K is available use the Epnp solver
    typedef openMVG::euclidean_resection::P3PSolver SolverType; // If K is available use the P3P solver
    MINIMUM_SAMPLES = SolverType::MINIMUM_SAMPLES;

    typedef ACKernelAdaptorResection_K<
      SolverType,  ResectionSquaredResidualError,  UnnormalizerResection, Mat34>  KernelType;

    KernelType kernel(pt2D, pt3D, *K);
    // Robustly estimation of the Projection matrix and it's precision
    std::pair<double,double> ACRansacOut = ACRANSAC(kernel, *pvec_inliers,
      ACRANSAC_ITER, P, dPrecision, true);
    *maxError = ACRansacOut.first;
  }

  // Test if the found model is valid
  if (pvec_inliers->size() > 2.5 * MINIMUM_SAMPLES)
  {
    // Re-estimate the model from the inlier data
    // and/or LM to Refine f R|t

    Mat3 K_, R_;
    Vec3 t_;
    KRt_From_P(*P, &K_, &R_, &t_);

    using namespace openMVG::bundle_adjustment;
    // Setup a BA problem
    BA_Problem_data<9> ba_problem;

    // Configure the size of the problem
    ba_problem.num_cameras_ = 1;
    ba_problem.num_points_ = pvec_inliers->size();
    ba_problem.num_observations_ = pvec_inliers->size();

    ba_problem.point_index_.reserve(ba_problem.num_observations_);
    ba_problem.camera_index_.reserve(ba_problem.num_observations_);
    ba_problem.observations_.reserve(2 * ba_problem.num_observations_);

    ba_problem.num_parameters_ = 9 * ba_problem.num_cameras_; // #[Rotation|translation|K] = [3x1]|[3x1]|[3x1]
    ba_problem.parameters_.reserve(ba_problem.num_parameters_);

    // Fill it with data (tracks and points coords)
    for (size_t i = 0; i < ba_problem.num_points_; ++i) {
      // Collect the image of point i in each frame.

      ba_problem.camera_index_.push_back(0);
      ba_problem.point_index_.push_back(i);
      const Vec2 & pt = pt2D.col((*pvec_inliers)[i]);
      ba_problem.observations_.push_back( pt(0) );
      ba_problem.observations_.push_back( pt(1) );
    }

    // Add camera parameters (R [3x1], t [3x1], focal [1])
    {
      // Rotation matrix to angle axis
      std::vector<double> angleAxis(3);
      ceres::RotationMatrixToAngleAxis((const double*) R_.data(), &angleAxis[0]);
      ba_problem.parameters_.push_back(angleAxis[0]);
      ba_problem.parameters_.push_back(angleAxis[1]);
      ba_problem.parameters_.push_back(angleAxis[2]);
      ba_problem.parameters_.push_back(t_[0]);
      ba_problem.parameters_.push_back(t_[1]);
      ba_problem.parameters_.push_back(t_[2]);
      ba_problem.parameters_.push_back(K_(0,0)); //focal
      ba_problem.parameters_.push_back(K_(0,2)); //principal_point x
      ba_problem.parameters_.push_back(K_(1,2)); //principal_point y
    }

    // Create residuals for each observation in the bundle adjustment problem. The
    // parameters for cameras and points are added automatically.
    ceres::Problem problem;
    for (size_t i = 0; i < ba_problem.num_points_; ++i) {
      // Each Residual block takes a point and a camera as input and outputs a 2
      // dimensional residual. Internally, the cost function stores the observed
      // image location and compares the reprojection against the observation.

      ceres::CostFunction* cost_function =
        new ceres::AutoDiffCostFunction<pinhole_reprojectionError::ErrorFunc_Refine_Camera, 2, 9>(
        new pinhole_reprojectionError::ErrorFunc_Refine_Camera(
        & ba_problem.observations()[2 * i],
        pt3D.col((*pvec_inliers)[i]).data()));

      problem.AddResidualBlock(cost_function,
        NULL, // squared loss
        ba_problem.mutable_camera_for_observation(0));
    }

    // Parameterization used to restrict camera intrinsics if desired (fixed focal length)
    ceres::SubsetParameterization *constant_transform_parameterization = NULL;
    if ( ! bRefineFocal ) {
      std::vector<int> vec_constant_focal(1,6); // Set last element camera parameter as fixed (the focal)
      constant_transform_parameterization =
        new ceres::SubsetParameterization(7, vec_constant_focal);
      problem.SetParameterization(ba_problem.mutable_camera_for_observation(0),
        constant_transform_parameterization);
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::SPARSE_SCHUR;
    if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::SUITE_SPARSE))
      options.sparse_linear_algebra_library_type = ceres::SUITE_SPARSE;
    else
      if (ceres::IsSparseLinearAlgebraLibraryTypeAvailable(ceres::CX_SPARSE))
        options.sparse_linear_algebra_library_type = ceres::CX_SPARSE;
      else
      {
        // No sparse back end for Ceres.
        // Use dense solving
        options.linear_solver_type = ceres::DENSE_SCHUR;
      }
    options.minimizer_progress_to_stdout = false;
    options.logging_type = ceres::SILENT;
#ifdef USE_OPENMP
    options.num_threads = omp_get_max_threads();
    options.num_linear_solver_threads = omp_get_max_threads();
#endif // USE_OPENMP

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // If no error, get back refined parameters

    if (summary.IsSolutionUsable())
    {
      const double * Rtf = ba_problem.mutable_camera_for_observation(0);

      // angle axis to rotation matrix
      ceres::AngleAxisToRotationMatrix(Rtf, R_.data());
      t_ = Vec3(Rtf[3], Rtf[4], Rtf[5]);
      double focal = Rtf[6];
      double ppx = Rtf[7];
      double ppy = Rtf[8];

      Mat3 KRefined;
      KRefined << focal,0, ppx,
        0, focal, ppy,
        0, 0, 1;

      PinholeCamera camRefined(KRefined,R_,t_);
      *P = camRefined._P;
    }

    return true;
  }
  else  {
    P = NULL;
    return false;
  }
}

} // namespace SfMRobust
} // namespace openMVG

#endif // OPENMVG_SFM_ROBUST_H

