
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_TRANSLATION_KNOWNROTATION_HPP
#define OPENMVG_MULTIVIEW_TRANSLATION_KNOWNROTATION_HPP

#include <vector>
#include "openMVG/multiview/projection.hpp"
#include "openMVG/multiview/solver_fundamental_kernel.hpp"
#include "openMVG/multiview/two_view_kernel.hpp"

//------------------
//-- Bibliography --
//------------------
//- [1] "Finding the Exact Rotation Between Two Images Independently of the Translation."
//- Authors: L.Kneip and R.Siegwart and M.Pollefeys.
//- Date: October 2012.
//- Conference: ECCV.

namespace openMVG {
namespace translation {
namespace kernel {

/**
 * Two-point translation estimation between two views from a known rotation
 * Implementation based on [1] => 3.2 Selection of the right solution. */
template<typename EpipolarDistanceErrorFunctor>
struct TwoPointTranslationSolver {
  enum { MINIMUM_SAMPLES = 2 };
  enum { MAX_MODELS = 1 };

  // Solve the problem of camera translation.
  static void Solve(const Mat &xA, const Mat &xB, const Mat3& R, std::vector<Vec3> *vec_t)
  {
    const Mat3 Rt = R.transpose();
    // A side bearing vectors
    const Vec3 f1 = Vec3(xA.col(0)(0), xA.col(0)(1), 1.);
    const Vec3 f2 = Vec3(xA.col(1)(0), xA.col(1)(1), 1.);
    // B side bearing vectors
    const Vec3 f1prime = Rt * Vec3(xB.col(0)(0), xB.col(0)(1), 1.);
    const Vec3 f2prime = Rt * Vec3(xB.col(1)(0), xB.col(1)(1), 1.);

    // Compute the translation of the camera
    const Vec3 c =
      ((f1.cross(f1prime))
      .cross(
      f2.cross(f2prime))).normalized();

    // Ensure the translation make the points in front of the cameras
    const Vec3 opticalFlow = f1 - f1prime;
    Vec3 translation = c;
    if (opticalFlow.dot(translation) < 0)
      translation = -translation;
    vec_t->push_back(-R * translation);
  }

  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  static double Error(const Mat3 &F, const Vec2 &x1, const Vec2 &x2){
  return EpipolarDistanceErrorFunctor::Error(F, x1, x2);
  }
};

//-- Generic Solver to find the translation from a known Rotation.
template<typename SolverArg,
  typename ErrorArg,
  typename ModelArg = Vec3>
class TranslationFromKnowRotation :
   public two_view::kernel::Kernel<SolverArg,ErrorArg, ModelArg>
{
public:
  // 2D / 2D points [camera coordinates]
  TranslationFromKnowRotation(const Mat &pt2DA, const Mat &pt2DB, const Mat3 & R) :
    two_view::kernel::Kernel<SolverArg, ErrorArg, ModelArg>(pt2DA, pt2DB), _R(R){}

  void Fit(const std::vector<size_t> &samples, std::vector<ModelArg> *models) const {
    const Mat ptA = ExtractColumns(this->x1_, samples);
    const Mat ptB = ExtractColumns(this->x2_, samples);

    assert(2 == ptA.rows());
    assert(2 == ptB.rows());

    SolverArg::Solve(ptA, ptB, _R, models);
  }

  // Error: distance of the sample to the epipolar line
  double Error(size_t sample, const ModelArg &model) const {
    Mat34 poseA, poseB;
    P_From_KRt(Mat3::Identity(), Mat3::Identity(), Vec3::Zero(), &poseA);
    P_From_KRt(Mat3::Identity(), _R, model, &poseB);
    const Mat3 F = F_from_P(poseA, poseB);
    return ErrorArg::Error(F, this->x1_.col(sample), this->x2_.col(sample));
  }
protected:
  const Mat3 _R;
};

//-- Usable solver for the 2pt translation from known rotation estimation
typedef TranslationFromKnowRotation<
  TwoPointTranslationSolver<openMVG::fundamental::kernel::SampsonError>, // SolverFunctor
  TwoPointTranslationSolver<openMVG::fundamental::kernel::SampsonError>, // ErrorFunctor
  Vec3>  TranslationFromKnowRotationKernel;

}  // namespace kernel
}  // namespace translation
}  // namespace openMVG

#endif  // OPENMVG_MULTIVIEW_TRANSLATION_KNOWNROTATION_HPP
