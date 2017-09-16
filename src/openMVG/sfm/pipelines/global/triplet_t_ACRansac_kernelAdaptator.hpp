// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013, 2014 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_GLOBAL_SFM_ENGINE_TRIPLET_T_ESTIMATOR_HPP
#define OPENMVG_GLOBAL_SFM_ENGINE_TRIPLET_T_ESTIMATOR_HPP

#include <vector>

#include "openMVG/multiview/conditioning.hpp"
#include "openMVG/numeric/numeric.h"

namespace openMVG{
namespace sfm{

/// AContrario Kernel to solve a translation triplet & structure problem
template <typename SolverArg,
          typename ErrorArg,
          typename ModelArg>
class TranslationTripletKernel_ACRansac
{
public:
  using Solver = SolverArg;
  using Model = ModelArg;

  TranslationTripletKernel_ACRansac
  (
    const Mat & x1,
    const Mat & x2,
    const Mat & x3,
    const std::vector<Mat3> & vec_KRi,
    const Mat3 & K,
    const double ThresholdUpperBound
  )
  : x1_(x1), x2_(x2), x3_(x3),
    Kinv_(K.inverse()),
    K_(K),
    logalpha0_(log10(M_PI)),
    ThresholdUpperBound_(ThresholdUpperBound),
    vec_KR_(vec_KRi)
  {
    // Normalize points by inverse(K)
    ApplyTransformationToPoints(x1_, Kinv_, &x1n_);
    ApplyTransformationToPoints(x2_, Kinv_, &x2n_);
    ApplyTransformationToPoints(x3_, Kinv_, &x3n_);

    vec_KR_[0] *= Kinv_;
    vec_KR_[1] *= Kinv_;
    vec_KR_[2] *= Kinv_;
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {

    // Create a model from the points
    Solver::Solve(ExtractColumns(x1n_, samples),
                  ExtractColumns(x2n_, samples),
                  ExtractColumns(x3n_, samples),
                  vec_KR_, models, ThresholdUpperBound_);
  }

  double Error(uint32_t sample, const Model &model) const {
    return ErrorArg::Error(model, x1n_.col(sample), x2n_.col(sample), x3n_.col(sample));
  }

  void Errors(const Model &model, std::vector<double> & vec_errors) const {
    for (size_t sample = 0; sample < x1n_.cols(); ++sample)
      vec_errors[sample] = ErrorArg::Error(model, x1n_.col(sample), x2n_.col(sample), x3n_.col(sample));
  }

  size_t NumSamples() const {
    return x1n_.cols();
  }

  void Unnormalize(Model * model) const {
    // Unnormalize model from the computed conditioning.
    model->P1 = K_ * model->P1;
    model->P2 = K_ * model->P2;
    model->P3 = K_ * model->P3;
  }

  double logalpha0() const {return logalpha0_;}

  double multError() const {return 1.0;}

  Mat3 normalizer1() const {return Kinv_;}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const { return sqrt(val) / Kinv_(0,0);}

private:
  const Mat & x1_, & x2_, & x3_;
  Mat x1n_, x2n_, x3n_;
  const Mat3 Kinv_, K_;
  const double logalpha0_;
  const double ThresholdUpperBound_;
  std::vector<Mat3> vec_KR_;

};

} // namespace sfm
} // namespace openMVG

#endif // OPENMVG_GLOBAL_SFM_ENGINE_TRIPLET_T_ESTIMATOR_HPP
