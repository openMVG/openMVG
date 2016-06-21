// Copyright (c) 2016 Pierre MOULON, Stephane FLOTRON

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "openMVG/linearProgramming/lInfinityCV/tijsAndXis_FromXi_Ri_non_central_cameras.hpp"
#include "openMVG/linearProgramming/linearProgramming.hpp"

namespace openMVG {
namespace sfm {
namespace non_central_camera {

/// AContrario Kernel to solve the non central cameras translation & structure problem
template
<
  typename SolverArg,
  typename ErrorArg,
  typename ModelArg>
class TranslationTripletKernel_ACRansac
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;

  TranslationTripletKernel_ACRansac
  (
    const std::vector< std::vector < std::vector <double> > > & pt,
    const std::vector<Mat3> & vec_global_rotations,
    const std::vector<Mat3> & vec_local_rotations,
    const std::vector<Vec3> & vec_local_centers,
    const double ThresholdUpperBound,
    const std::pair<IndexT, IndexT> image_dimension
  )
  : pt_(pt),
    vec_global_rotations_(vec_global_rotations),
    vec_local_rotations_(vec_local_rotations),
    vec_local_centers_(vec_local_centers),
    ThresholdUpperBound_(ThresholdUpperBound),
    logalpha0_(log10(M_PI/(image_dimension.first*image_dimension.second)))
  {
    // loagalpha configured to use pixel residual error units
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit
  (
    const std::vector<size_t> &samples,
    std::vector<Model> *models
  ) const
  {
    std::vector < std::vector < std::vector <double > > > pt_sampled(samples.size());

    for (int i=0; i < samples.size(); ++i)
      pt_sampled[i] = pt_[samples[i]];

    // Create a model from the points
    Solver::Solve(
      pt_sampled,
      vec_local_rotations_, vec_local_centers_,
      vec_global_rotations_, models, ThresholdUpperBound_);
  }

  void Errors
  (
    const Model &model,
    std::vector<double> & vec_errors
  ) const
  {
    for (size_t sample = 0; sample < pt_.size(); ++sample)
      vec_errors[sample] = ErrorArg::Error(model, pt_[sample], vec_local_rotations_, vec_local_centers_);
  }

  double Error
  (
    size_t sample,
    const Model &model
  ) const
  {
    return ErrorArg::Error(model, pt_[sample], vec_local_rotations_, vec_local_centers_);
  }

  size_t NumSamples() const
  {
    return pt_.size();
  }

  void Unnormalize(Model * model) const
  {
    // Unnormalize model is not necessary since normalized point are used
  }

  Mat3 normalizer1() const {return Mat3::Identity();}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const { return val;}

  double logalpha0() const {return logalpha0_;}

  double multError() const {return 1.0;}

private:
  const std::vector < std::vector < std::vector <double > > > & pt_;
  const double logalpha0_;
  const double ThresholdUpperBound_;
  const std::vector<Mat3> & vec_global_rotations_;
  const std::vector<Mat3> & vec_local_rotations_;
  const std::vector<Vec3> & vec_local_centers_;

};
} // namespace non_central_camera
} // namespace sfm
} // namespace openMVG
