#pragma once

#include <openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp>
#include <openMVG/numeric/numeric.h>
#include <openMVG/multiview/conditioning.hpp>

namespace openMVG {
namespace robust {


template <typename SolverArg,
  typename ErrorArg,
  typename UnnormalizerArg,
  typename ModelArg = Mat3,
  typename SolverLSArg = SolverArg>
class KernelAdaptorLoRansac : public ACKernelAdaptor<SolverArg, ErrorArg, UnnormalizerArg, ModelArg>
{
public:
  typedef SolverArg Solver;
  typedef ModelArg Model;
  typedef ErrorArg ErrorT;
  typedef SolverLSArg SolverLS;
  
  enum
  {
    MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES,
    MINIMUM_LSSAMPLES = SolverLS::MINIMUM_SAMPLES
  };

  KernelAdaptorLoRansac(const Mat &x1, int w1, int h1,
                        const Mat &x2, int w2, int h2, 
                        bool bPointToLine = true) 
          : ACKernelAdaptor<SolverArg, ErrorArg, UnnormalizerArg, ModelArg>(x1, w1, h1, x2, w2, h2, bPointToLine)
  {
  }
  
  void FitLS(const std::vector<std::size_t> &inliers, 
              std::vector<Model> *models, 
              const std::vector<double> *weights = nullptr) const
  {
    const Mat x1 = ExtractColumns(this->x1_, inliers);
    const Mat x2 = ExtractColumns(this->x2_, inliers);
    SolverLS::Solve(x1, x2, models, weights);
  }
  
  void computeWeights(const Model & model, 
                      const std::vector<std::size_t> &inliers, 
                      std::vector<double> & vec_weights, 
                      const double eps = 0.001) const
  {
    const auto numInliers = inliers.size();
    vec_weights.resize(numInliers);
    for(std::size_t sample = 0; sample < numInliers; ++sample)
    {
      const auto idx = inliers[sample];
      vec_weights[sample] = ErrorT::Error(model, this->x1_.col(idx), this->x2_.col(idx));
      // avoid division by zero
      vec_weights[sample] = 1.0 / std::pow(std::max(eps, vec_weights[sample]), 2.0);
    }
  }

};


template <typename SolverArg,
typename ErrorArg,
typename UnnormalizerArg,
typename SolverLSArg,
typename ModelArg = Mat34>
class KernelAdaptorResectionLORansac_K : 
      public ACKernelAdaptorResection_K<SolverArg, ErrorArg, UnnormalizerArg, ModelArg>
{
  public:
  typedef SolverArg Solver;
  typedef ModelArg Model;
  typedef ErrorArg ErrorT;
  typedef SolverLSArg SolverLS;
  
  enum
  {
    MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES,
    MINIMUM_LSSAMPLES = SolverLS::MINIMUM_SAMPLES
  };
  
  KernelAdaptorResectionLORansac_K(const Mat &x2d, const Mat &x3D, const Mat3 & K) : 
      ACKernelAdaptorResection_K<SolverArg, ErrorArg, UnnormalizerArg, ModelArg>(x2d, x3D, K)
  {}
  
  void FitLS(const std::vector<std::size_t> &inliers, 
              std::vector<Model> *models, 
              const std::vector<double> *weights = nullptr) const
  {
    const Mat x1 = ExtractColumns(this->x2d_, inliers);
    const Mat x2 = ExtractColumns(this->x3D_, inliers);
    SolverLS::Solve(x1, x2, models, weights);
  }
  
  void computeWeights(const Model & model, 
                      const std::vector<std::size_t> &inliers, 
                      std::vector<double> & vec_weights, 
                      const double eps = 0.001) const
  {
    const auto numInliers = inliers.size();
    vec_weights.resize(numInliers);
    for(std::size_t sample = 0; sample < numInliers; ++sample)
    {
      const auto idx = inliers[sample];
      vec_weights[sample] = ErrorT::Error(model, this->x2d_.col(idx), this->x3D_.col(idx));
      // avoid division by zero
      vec_weights[sample] = 1/std::pow(std::max(eps, vec_weights[sample]), 2);
    }
  }
};




} // namespace robust
} // namespace openMVG