// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.
//
//:\file
//\author Ricardo Fabbri, Brown & Rio de Janeiro State U. (rfabbri.github.io) 
//\date Tue Jun  1 11:55:58 -03 2021
//\author Gabriel ANDRADE Rio de Janeiro State U.
//\author Pierre MOULON
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_MULTIVIEW_THREE_VIEW_KERNEL_HPP
#define OPENMVG_MULTIVIEW_THREE_VIEW_KERNEL_HPP

#include <vector>
#include "openMVG/numeric/extract_columns.hpp"
#include "openMVG/multiview/multiview_match_constraint.hpp"
#include "openMVG/multiview/trifocal/trifocal_model.hpp"

namespace openMVG {
namespace trifocal {
  
template<typename SolverArg,
         typename ErrorArg,
         typename ModelArg = trifocal_model_t>
class ThreeViewKernel {
public:
  using Solver = SolverArg;
  using Model = ModelArg;
  using ErrorT = ErrorArg;
  /// The minimal number of point required for the model estimation
  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  /// The number of models that the minimal solver could return.
  enum { MAX_MODELS = Solver::MAX_MODELS };
  
  ThreeViewKernel(const Mat &x1, const Mat &x2, const Mat &x3) 
    : x1_(x1), x2_(x2), x3_(x3), multiview_match_constraint_(MultiviewMatchConstraint::ORIENTATION) {}

  /// Extract required sample and fit model(s) to the sample
  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {
    const auto
      x1 = ExtractColumns(x1_, samples),
      x2 = ExtractColumns(x2_, samples),
      x3 = ExtractColumns(x3_, samples);
    Solver::Solve(x1, x2, x3, models);
  }
  
  /// Return the error associated to the model and sample^nth point
  double Error(uint32_t sample, const Model &model) const {
    return ErrorArg::Error(model, x1_.col(sample), x2_.col(sample), x3_.col(sample), multiview_match_constraint_);
  }

  /// Number of putative point
  size_t NumSamples() const { return static_cast<size_t>(x1_.cols()); }

  /// Compute a model on sampled datum_
  static void Solve(const Mat &x1, const Mat &x2, const Mat &x3, std::vector<Model> *models) {
    Solver::Solve(x1, x2, x3, models); // By offering this, Kernel types can be passed to templates.
  }

  void SetMultiviewMatchConstraint(MultiviewMatchConstraint c) 
  { multiview_match_constraint_ = c; }
  
protected:
  const Mat &x1_, &x2_, &x3_; // corresponding point of the trifocal configuration
  MultiviewMatchConstraint multiview_match_constraint_;
};

} // namespace trifocal
} // namespace OpenMVG
#endif  // OPENMVG_MULTIVIEW_THREE_VIEW_KERNEL_HPP
