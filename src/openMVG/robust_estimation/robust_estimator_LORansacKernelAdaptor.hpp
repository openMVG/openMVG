#pragma once

#include <openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp>
#include <openMVG/numeric/numeric.h>
#include <openMVG/multiview/conditioning.hpp>

namespace openMVG {
namespace robust {


/**
 * @brief The generic kernel to be used for LORansac framework.
 * 
 * @tparam SolverArg The minimal solver able to find a solution from a
 * minimum set of points.
 * @tparam ErrorArg The functor computing the error for each data sample with
 * respect to the estimated model.
 * @tparam UnnormalizerArg The functor used to normalize the data before the 
 * estimation of the model.
 * @tparam ModelArg = Mat3 The type of the model to estimate.
 * @tparam SolverLSArg = SolverArg The least square solver that is used to find
 * a solution from any set of data larger than the minimum required.
 */
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
  
  /**
   * @brief Given a model and the associated inliers, it computes the weight for
   * each inlier as the squared inverse of the associated error.
   * 
   * @param[in] model The model to against which to compute the weights.
   * @param[in] inliers The inliers associated to the model.
   * @param[out] vec_weights The weights associated to each inlier.
   * @param[in] eps Each inlier having an error below this value will be assigned
   * a weight of 1/eps^2 (to avoid division by zero).
   */
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

/**
 * @brief The kernel for the resection with known intrinsics (PnP) to be used with
 * the LORansac framework.
 * 
 * @tparam SolverArg The minimal solver able to find a solution from a
 * minimum set of points, usually any PnP solver.
 * @tparam ErrorArg The functor computing the error for each data sample with
 * respect to the estimated model, usually a reprojection error functor.
 * @tparam UnnormalizerArg The functor used to normalize the data before the 
 * estimation of the model, usually a functor that normalize the point in camera
 * coordinates (ie multiply by the inverse of the calibration matrix).
 * @tparam ModelArg = Mat34 The type of the model to estimate, the projection matrix.
 * @tparam SolverLSArg = SolverArg The least square solver that is used to find
 * a solution from any set of data larger than the minimum required, usually the 
 * 6 point algorithm which solves the resection problem by means of LS.
 */
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