// Copyright (c) 2016 Pierre MOULON, Stephane FLOTRON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <vector>
#include "openMVG/numeric/numeric.h"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"

#include "opengv/types.hpp"
#include "opengv/relative_pose/methods.hpp"
#include "opengv/triangulation/methods.hpp"
#include "opengv/relative_pose/NoncentralRelativeAdapter.hpp"
#include "opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp"
#include <opengv/sac/Ransac.hpp>

/*************************
 *
 *  OpenMVG/OpenGV Wrapper for relative pose computation between non central cameras
 *  - Kernel adaptor for camera rig relative poses computation:
 *    - Relative pose wrapper (GE/SIXPT)
 *    - Non central cameras pixel/angular reprojection error computation
 *  - A Contrario Kernel adaptor for the non central
 *  - Kernel adaptor for camera rig relative poses computation
 *
 ************************
 */

namespace openMVG{
namespace noncentral{
namespace kernel{

using namespace opengv;

/**
 * Six point solver for non central camera system,
 * // [1] "Solutions to minimal generalized relative pose problems".
 * // authors: Stewenius, H., Nister, D., Oskarsson, M., & Astrom K,
 * // Date: 2005:
 * // Conference: Workshop on omnidirectional vision 2005.
 */
struct SixPointSolver
{
  enum { MINIMUM_SAMPLES = 9 }; // /!\ Using more point provide better results
  enum { MAX_MODELS = 1 };
  static void Solve
  (
    relative_pose::NoncentralRelativeAdapter & adapter,
    std::vector<transformation_t> * models,
    const std::vector<size_t> &indices
  )
  {
    // convert size_t to int for opengv call
    const std::vector<int> idx(indices.begin(), indices.end());

    // create non central relative sac problem
    sac_problems::relative_pose::NoncentralRelativePoseSacProblem problem(
      adapter,
      sac_problems::relative_pose::NoncentralRelativePoseSacProblem::SIXPT,
      false);

    // solve pose problem
    transformation_t relativePose;
    if (problem.computeModelCoefficients(idx, relativePose))
    {
      // Refine to compute the translation
      adapter.sett12(relativePose.col(3));
      adapter.setR12(relativePose.block<3,3>(0,0));
      relativePose = relative_pose::optimize_nonlinear(adapter);

      models->push_back(relativePose);
    }
  }
};

/**
 * Generalized eigen value solver for non central camera system,
 * @InProceedings{ Kneip_2014_CVPR,
 * author = {Kneip, Laurent and Li, Hongdong},
 * title = {Efficient Computation of Relative Pose for Multi-Camera Systems},
 * journal = {The IEEE Conference on Computer Vision and Pattern Recognition (CVPR)},
 * month = {June},
 * year = {2014}
 * }
 */

struct GePointSolver
{
  enum { MINIMUM_SAMPLES = 6 };
  enum { MAX_MODELS = 1 };
  static void Solve
  (
    relative_pose::NoncentralRelativeAdapter & adapter,
    std::vector<transformation_t> * models,
    const std::vector<size_t> &indices
  )
  {
    // convert size_t to int for opengv call
    const std::vector<int> idx(indices.begin(), indices.end());

    // create non central relative sac problem
    sac_problems::relative_pose::NoncentralRelativePoseSacProblem problem(
      adapter,
      sac_problems::relative_pose::NoncentralRelativePoseSacProblem::GE,
      false);

    // solve pose problem
    transformation_t relativePose;
    if (problem.computeModelCoefficients(idx, relativePose))
    {
      // Refine the model
      adapter.sett12(relativePose.col(3));
      adapter.setR12(relativePose.block<3,3>(0,0));
      relativePose = relative_pose::optimize_nonlinear(adapter, idx);

      models->push_back(relativePose);
    }
  }
};

// compute reprojection error
struct RigProjError
{
  static double Error
  (
    size_t sample,
    const transformation_t & relativePose,
    relative_pose::NoncentralRelativeAdapter & adapter
  )
  {
    // extract pose of cameras
    const Vec3 bearingOne = adapter.getBearingVector1(sample);
    const Vec3 bearingTwo = adapter.getBearingVector2(sample);

    const Vec2 x1 = bearingOne.head(2) / bearingOne(2);
    const Vec2 x2 = bearingTwo.head(2) / bearingTwo(2);

    const Mat3 R1 = adapter.getCamRotation1(sample).transpose();
    const Mat3 R2 = adapter.getCamRotation2(sample).transpose();

    const Vec3 t1 = - R1 * adapter.getCamOffset1(sample);
    const Vec3 t2 = - R2 * adapter.getCamOffset2(sample);

    // retrieve relative pose of rigs
    const translation_t CRig = relativePose.col(3);
    const rotation_t rRig = relativePose.block<3,3>(0,0).transpose();
    const Vec3  tRig = -rRig * CRig;

    // compute relative pose of cameras
    const rotation_t R = R2 * rRig ;
    const translation_t t = R2 * tRig + t2 ;

    // compute 3d point and reprojection error
    const Mat34 P1 = HStack(R1, t1);
    const Mat34 P2 = HStack(R, t);

    // Triangulate and return the reprojection error
    Triangulation triangulationObj;
    triangulationObj.add(P1, x1);
    triangulationObj.add(P2, x2);

    const Vec3 X = triangulationObj.compute();

    //- Return max error as a test
    const double pt1ReProj = (Project(P1, X) - x1).norm();
    const double pt2ReProj = (Project(P2, X) - x2).norm();

    return std::max(pt1ReProj, pt2ReProj);
  }
};

// compute angular error (as in kneip opengv library)
struct RigAngularError
{
  static double Error
  (
    size_t index,
    const transformation_t & model,
    relative_pose::NoncentralRelativeAdapter & adapter
  )
  {
    // extract rotation and translation from model
    const translation_t translation = model.col(3);
    const rotation_t rotation = model.block<3,3>(0,0);

    // compute pose
    const translation_t cam1Offset = adapter.getCamOffset1(index);
    const rotation_t cam1Rotation = adapter.getCamRotation1(index);
    const translation_t cam2Offset = adapter.getCamOffset2(index);
    const rotation_t cam2Rotation = adapter.getCamRotation2(index);

    const translation_t directTranslation =
        cam1Rotation.transpose() *
        ((translation - cam1Offset) + rotation * cam2Offset);
    const rotation_t directRotation =
        cam1Rotation.transpose() * rotation * cam2Rotation;

    adapter.sett12(directTranslation);
    adapter.setR12(directRotation);

    transformation_t inverseSolution;
    inverseSolution.block<3,3>(0,0) = directRotation.transpose();
    inverseSolution.col(3) =
        -inverseSolution.block<3,3>(0,0)*directTranslation;

    Vec4 p_hom;
    p_hom << opengv::triangulation::triangulate2(adapter,index), 1.0;

    const bearingVector_t reprojection1 = p_hom.block<3,1>(0,0).normalized();
    const bearingVector_t reprojection2 = (inverseSolution * p_hom).normalized();
    const bearingVector_t f1 = adapter.getBearingVector1(index).normalized();
    const bearingVector_t f2 = adapter.getBearingVector2(index).normalized();

    //bearing-vector based outlier criterium (select threshold accordingly):
    //1-(f1'*f2) = 1-cos(alpha) \in [0:2]
    const double reprojError1 = 1.0 - f1.transpose() * reprojection1;
    const double reprojError2 = 1.0 - f2.transpose() * reprojection2;

    return std::max(reprojError1,reprojError2);
  }
};

}  // namespace kernel

namespace robust{

 using namespace opengv;
  /// Rig Pose Kernel adaptator for the A contrario model estimator
  template
    <
      typename SolverArg,
      typename ErrorArg,
      typename ModelArg = transformation_t
    >
  class ACKernelAdaptor_non_central_relative_pose
  {
  public:
    typedef SolverArg Solver;
    typedef ModelArg  Model;
    typedef ErrorArg ErrorT;

    ACKernelAdaptor_non_central_relative_pose
    (
      const bearingVectors_t & b1,
      const bearingVectors_t & b2,
      const std::vector<int> & scIdOne,
      const std::vector<int> & scIdTwo,
      const translations_t & rigOffsets,
      const rotations_t & rigRotations
    )
    : b1_(b1), b2_(b2),
      scIdOne_(scIdOne), scIdTwo_(scIdTwo),
      offSets(rigOffsets), rotations(rigRotations),
      logalpha0_(log10(0.5))
    {
        assert(b1_.size() == b2_.size());
    }

    enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
    enum { MAX_MODELS = Solver::MAX_MODELS };

    void Fit
    (
      const std::vector<size_t> &samples,
      std::vector<Model> *models
    )
    const
    {
      //create non-central relative adapter
      relative_pose::NoncentralRelativeAdapter adapter(
            b1_, b2_, scIdOne_ , scIdTwo_,
            offSets, rotations);

      Solver::Solve(adapter, models, samples);
    }

    void Errors
    (
      const Model & model,
      std::vector<double> & vec_errors
    )
    const
    {
      //create non-central relative adapter
      relative_pose::NoncentralRelativeAdapter adapter(
            b1_, b2_, scIdOne_ , scIdTwo_,
            offSets, rotations);

      vec_errors.resize(b1_.size());
      for (size_t sample = 0; sample < b1_.size(); ++sample)
        vec_errors[sample] = ErrorT::Error(sample, model, adapter);
    }

    double Error
    (
      size_t sample,
      const Model & model
    )
    const
    {
      //create non-central relative adapter
      relative_pose::NoncentralRelativeAdapter adapter(
            b1_, b2_, scIdOne_ , scIdTwo_,
            offSets, rotations);

      return ErrorT::Error(sample, model, adapter);
    }

    size_t NumSamples() const { return b1_.size(); }
    void Unnormalize(Model * model) const {}
    double logalpha0() const {return logalpha0_;}
    double multError() const {return 0.5;} // angle vs. angle error
    Mat3 normalizer1() const {return Mat3::Identity();}
    Mat3 normalizer2() const {return Mat3::Identity();}
    double unormalizeError(double val) const { return val; }

  private:
    const bearingVectors_t & b1_, b2_; // bearing vectors.
    const std::vector<int> & scIdOne_, scIdTwo_ ; // cam correspondences
    const translations_t & offSets ; // camera position in rigid rig frame
    const rotations_t & rotations ; // rig subcamera orientation
    double logalpha0_; // Alpha0 is used to make the error adaptive to the image size
  };
} // end robust

namespace robust_relative_pose{

using namespace openMVG::matching;

/**
 * @brief Estimate the relative pose between two non central camera from corresponding bearing vectors
 *
 * @param[in] b1 bearing vectors of the non central camera 1
 * @param[in] b2 bearing vectors of the non central camera 2
 * @param[in] scIdOne subcamera id of each bearing vector of the non central camera 1
 * @param[in] scIdTwo subcamera id of each bearing vector of the non central camera 2
 * @param[in] rigOffsets non central camera centers (local coordinates)
 * @param[in] rigRotations non central camera rotations (local coordinates)
 * @param[out] transformation_t 'found' relative pose (R^T and C)
 * @param[out] pvec_inliers 'found 'inliers indices
 * @param[out] errorMax upper bound of the reprojection error of the found solution
 * @param[in] precision upper bound of the desired solution
 */
bool non_central_cameras_robust_relative_pose
(
  const opengv::bearingVectors_t & b1,
  const opengv::bearingVectors_t & b2,
  const std::vector<int> & scIdOne,
  const std::vector<int> & scIdTwo,
  const opengv::translations_t & rigOffsets,
  const opengv::rotations_t & rigRotations,
  opengv::transformation_t * relativePose,
  std::vector<size_t> * pvec_inliers,
  double * errorMax,
  double precision
)
{
  assert(pvec_inliers != NULL);

  // Use the Generalized eigenvalue solver to solve pose problem
  typedef openMVG::noncentral::kernel::GePointSolver SolverType;

  // Define the AContrario adaptor
  typedef noncentral::robust::ACKernelAdaptor_non_central_relative_pose
    <
      SolverType,
      openMVG::noncentral::kernel::RigAngularError,
      opengv::transformation_t
    > KernelType;

  KernelType kernel(b1, b2, scIdOne, scIdTwo, rigOffsets, rigRotations);

  // Robustly estimation of the Relative pose matrix and it's associated AContrario precision
  const std::pair<double,double> acRansacOut =
    openMVG::robust::ACRANSAC(
      kernel,
      *pvec_inliers,
      4096,
      relativePose,
      precision,
      false);

  if (errorMax)
    *errorMax = acRansacOut.first;

  return
    (
       pvec_inliers->size() > 2.5 * SolverType::MINIMUM_SAMPLES * rigOffsets.size() // Sufficient inlier coverage
       && acRansacOut.first != std::numeric_limits<double>::infinity() // Check AContrario estimation found a model
    );

  // If you want to use the OpenGV RANSAC, please consider this version
  // Use it will be slower than using AC-RANSAC

//  using namespace opengv;
//  relative_pose::NoncentralRelativeAdapter adapter(
//    b1, b2, scIdOne, scIdTwo, rigOffsets, rigRotations);
//
//  sac::Ransac<
//      sac_problems::relative_pose::NoncentralRelativePoseSacProblem> ransac;
//  std::shared_ptr<
//      sac_problems::relative_pose::NoncentralRelativePoseSacProblem>
//      relposeproblem_ptr(
//      new sac_problems::relative_pose::NoncentralRelativePoseSacProblem(
//      adapter,
//      sac_problems::relative_pose::NoncentralRelativePoseSacProblem::GE));
//  ransac.sac_model_ = relposeproblem_ptr;
//  ransac.threshold_ = precision;
//  ransac.max_iterations_ = 4096;
//
//  ransac.computeModel();
//
//  *relativePose = ransac.model_coefficients_;
//
//  return pvec_inliers->size() > 2.5 * 6 * rigOffsets.size();
}

} // robust_relative_pose
} // namespace noncentral
} // namespace openMVG
