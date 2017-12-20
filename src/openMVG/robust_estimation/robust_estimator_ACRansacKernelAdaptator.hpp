// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_KERNEL_ADAPTATOR_HPP
#define OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_KERNEL_ADAPTATOR_HPP

// Here a collection of A contrario Kernel adaptor.
//  - See // [1] "Robust and accurate calibration of camera networks". PhD.
//  - Authors: Pierre MOULON
//
//ACKernelAdaptor
//  i.e. general two view estimation (affine, homography, fundamental)
//
// For pose/resection estimation:
//  - ACKernelAdaptorResection
// For essential matrix estimation
//  - ACKernelAdaptorEssential
//  - ACKernelAdaptor_AngularRadianError
//
// Mainly it add correct data normalization and define the required functions
//  by the ACRANSAC algorithm.
//

#include <vector>

#include "openMVG/multiview/conditioning.hpp"
#include "openMVG/numeric/extract_columns.hpp"

namespace openMVG {
namespace robust{

/// Two view Kernel adapter for the A contrario model estimator
/// Handle data normalization and compute the corresponding logalpha 0
///  that depends of the error model (point to line, or point to point)
/// This kernel adapter is working for affine, homography, fundamental matrix
///  estimation.
template <typename SolverArg,
          typename ErrorArg,
          typename UnnormalizerArg,
          typename ModelArg = Mat3>
class ACKernelAdaptor
{
public:
  using Solver = SolverArg;
  using Model = ModelArg;
  using ErrorT = ErrorArg;

  ACKernelAdaptor(
    const Mat &x1, int w1, int h1,
    const Mat &x2, int w2, int h2, bool bPointToLine = true)
    : x1_(x1.rows(), x1.cols()), x2_(x2.rows(), x2.cols()),
    N1_(3,3), N2_(3,3), logalpha0_(0.0), bPointToLine_(bPointToLine)
  {
    assert(2 == x1_.rows());
    assert(x1_.rows() == x2_.rows());
    assert(x1_.cols() == x2_.cols());

    NormalizePoints(x1, &x1_, &N1_, w1, h1);
    NormalizePoints(x2, &x2_, &N2_, w2, h2);

    // LogAlpha0 is used to make error data scale invariant
    if (bPointToLine)  {
      // Ratio of containing diagonal image rectangle over image area
      const double D = std::hypot(w2, h2); // diameter
      const double A = w2*(double)h2; // area
      logalpha0_ = log10(2.0*D/A /N2_(0,0));
    }
    else  {
      // ratio of area : unit circle over image area
      logalpha0_ = log10(M_PI/(w2*(double)h2) /(N2_(0,0)*N2_(0,0)));
    }
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {
    const Mat
      x1 = ExtractColumns(x1_, samples),
      x2 = ExtractColumns(x2_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(uint32_t sample, const Model &model) const {
    return ErrorT::Error(model, x1_.col(sample), x2_.col(sample));
  }

  void Errors(const Model & model, std::vector<double> & vec_errors) const
  {
    vec_errors.resize(x1_.cols());
    for (uint32_t sample = 0; sample < x1_.cols(); ++sample)
      vec_errors[sample] = ErrorT::Error(model, x1_.col(sample), x2_.col(sample));
  }

  size_t NumSamples() const {
    return static_cast<size_t>(x1_.cols());
  }

  void Unnormalize(Model * model) const {
    // Unnormalize model from the computed conditioning.
    UnnormalizerArg::Unnormalize(N1_, N2_, model);
  }

  double logalpha0() const {return logalpha0_;}

  double multError() const {return (bPointToLine_)? 0.5 : 1.0;}

  Mat3 normalizer1() const {return N1_;}
  Mat3 normalizer2() const {return N2_;}
  double unormalizeError(double val) const {return sqrt(val) / N2_(0,0);}

private:
  Mat x1_, x2_;       // Normalized input data
  Mat3 N1_, N2_;      // Matrix used to normalize data
  double logalpha0_;  // Alpha0 is used to make the error adaptive to the image size
  bool bPointToLine_; // Store if error model is pointToLine or point to point
};

struct UnnormalizerResection {
  static void Unnormalize(const Mat &T, const Mat &U, Mat34 *P){
    *P = T.inverse() * (*P);
  }
};

/// Pose/Resection Kernel adapter for the A contrario model estimator
template <typename SolverArg,
  typename ErrorArg,
  typename UnnormalizerArg,
  typename ModelArg = Mat34>
class ACKernelAdaptorResection
{
public:
  using Solver = SolverArg;
  using Model = ModelArg;
  using ErrorT = ErrorArg;

  ACKernelAdaptorResection(const Mat &x2d, int w, int h, const Mat &x3D)
    : x2d_(x2d.rows(), x2d.cols()), x3D_(x3D),
    N1_(3,3), logalpha0_(log10(M_PI))
  {
    assert(2 == x2d_.rows());
    assert(3 == x3D_.rows());
    assert(x2d_.cols() == x3D_.cols());

    NormalizePoints(x2d, &x2d_, &N1_, w, h);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {
    const Mat
      x1 = ExtractColumns(x2d_, samples),
      x2 = ExtractColumns(x3D_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(uint32_t sample, const Model &model) const {
    return ErrorT::Error(model, x2d_.col(sample), x3D_.col(sample));
  }

  void Errors(const Model & model, std::vector<double> & vec_errors) const
  {
    vec_errors.resize(x2d_.cols());
    for (uint32_t sample = 0; sample < x2d_.cols(); ++sample)
      vec_errors[sample] = ErrorT::Error(model, x2d_.col(sample), x3D_.col(sample));
  }

  size_t NumSamples() const { return x2d_.cols(); }

  void Unnormalize(Model * model) const {
    // Unnormalize model from the computed conditioning.
    UnnormalizerArg::Unnormalize(N1_, Mat3::Identity(), model);
  }

  double logalpha0() const {return logalpha0_;}
  double multError() const {return 1.0;} // point to point error
  Mat3 normalizer1() const {return Mat3::Identity();}
  Mat3 normalizer2() const {return N1_;}
  double unormalizeError(double val) const {return sqrt(val) / N1_(0,0);}

private:
  Mat x2d_;
  const Mat & x3D_;
  Mat3 N1_;      // Matrix used to normalize data
  double logalpha0_; // Alpha0 is used to make the error adaptive to the image size
};

/// Essential matrix Kernel adaptor for the A contrario model estimator
template <typename SolverArg,
  typename ErrorArg,
  typename ModelArg = Mat3>
class ACKernelAdaptorEssential
{
public:
  using Solver = SolverArg;
  using Model = ModelArg;
  using ErrorT = ErrorArg;

  ACKernelAdaptorEssential
  (
    const Mat2X &x1, const Mat3X & bearing1, int w1, int h1,
    const Mat2X &x2, const Mat3X & bearing2, int w2, int h2,
    const Mat3 & K1, const Mat3 & K2
  ):x1_(x1),
    x2_(x2),
    bearing1_(bearing1),
    bearing2_(bearing2),
    N1_(Mat3::Identity()),
    N2_(Mat3::Identity()), logalpha0_(0.0),
    K1_(K1), K2_(K2)
  {
    assert(2 == x1_.rows());
    assert(x1_.rows() == x2_.rows());
    assert(x1_.cols() == x2_.cols());

    //Point to line probability (line is the epipolar line)
    const double D = std::hypot(w2, h2); // diameter
    const double A = w2*(double)h2; // area
    logalpha0_ = log10(2.0*D/A * .5);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {
    const Mat
      x1 = ExtractColumns(bearing1_, samples),
      x2 = ExtractColumns(bearing2_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(uint32_t sample, const Model &model) const {
    Mat3 F;
    FundamentalFromEssential(model, K1_, K2_, &F);
    return ErrorT::Error(F, this->x1_.col(sample), this->x2_.col(sample));
  }

  void Errors(const Model & model, std::vector<double> & vec_errors) const
  {
    Mat3 F;
    FundamentalFromEssential(model, K1_, K2_, &F);
    vec_errors.resize(x1_.cols());
    for (uint32_t sample = 0; sample < x1_.cols(); ++sample)
      vec_errors[sample] = ErrorT::Error(F, this->x1_.col(sample), this->x2_.col(sample));
  }

  size_t NumSamples() const { return x1_.cols(); }
  void Unnormalize(Model * model) const {}
  double logalpha0() const {return logalpha0_;}
  double multError() const {return 0.5;} // point to line error
  Mat3 normalizer1() const {return N1_;}
  Mat3 normalizer2() const {return N2_;}
  double unormalizeError(double val) const { return val; }

private:
  Mat2X x1_, x2_;       // image points
  Mat3X bearing1_, bearing2_;   // bearing vectors
  Mat3 N1_, N2_;      // Matrix used to normalize data
  double logalpha0_;  // Alpha0 is used to make the error adaptive to the image size
  Mat3 K1_, K2_;      // Intrinsic camera parameter
};

/// Two view Kernel adapter for the A contrario model estimator.
/// Specialization to handle radian angular residual error.
template <typename SolverArg,
          typename ErrorArg,
          typename ModelArg = Mat3>
class ACKernelAdaptor_AngularRadianError
{
public:
  using Solver = SolverArg;
  using Model = ModelArg;
  using ErrorT = ErrorArg;

  ACKernelAdaptor_AngularRadianError(
    const Mat & xA,
    const Mat & xB):
    x1_(xA), x2_(xB),
    logalpha0_(log10(1.0/2.0))
  {
    assert(3 == x1_.rows());
    assert(x1_.rows() == x2_.rows());
    assert(x1_.cols() == x2_.cols());
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {
    const Mat
      x1 = ExtractColumns(x1_, samples),
      x2 = ExtractColumns(x2_, samples);
    Solver::Solve(x1, x2, models);
  }

  double Error(uint32_t sample, const Model &model) const {
    return Square(ErrorT::Error(model, x1_.col(sample), x2_.col(sample)));
  }

  void Errors(const Model & model, std::vector<double> & vec_errors) const
  {
    vec_errors.resize(x1_.cols());
    for (uint32_t sample = 0; sample < x1_.cols(); ++sample)
      vec_errors[sample] = Square(ErrorT::Error(model, x1_.col(sample), x2_.col(sample)));
  }

  size_t NumSamples() const {
    return static_cast<size_t>(x1_.cols());
  }

  void Unnormalize(Model * model) const {
    //-- Do nothing, no normalization in the angular case
  }

  double logalpha0() const {return logalpha0_;}

  double multError() const {return 1./4.;}

  Mat3 normalizer1() const {return Mat3::Identity();}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const {return sqrt(val);}

private:
  Mat x1_, x2_;      // Normalized input data
  double logalpha0_; // Alpha0 is used to make the error scale invariant
};

} // namespace robust
} // namespace openMVG

#endif // OPENMVG_ROBUST_ESTIMATOR_ACRANSAC_KERNEL_ADAPTATOR_HPP
