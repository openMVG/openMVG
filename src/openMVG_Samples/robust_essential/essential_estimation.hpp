
#ifndef OPENMVG_SAMPLES_ROBUST_ESSENTIAL_HPP
#define OPENMVG_SAMPLES_ROBUST_ESSENTIAL_HPP

#include "openMVG/multiview/solver_essential_kernel.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansacKernelAdaptator.hpp"

using namespace openMVG::robust;

namespace openMVG
{

/// Estimate the essential matrix from point matches and K matrices.
bool robustEssential(
  const Mat3 & K1, const Mat3 & K2,
  const Mat & x1, const Mat & x2,
  Mat3 * pE,
  std::vector<size_t> * pvec_inliers,
  const std::pair<size_t, size_t> & size_ima1,
  const std::pair<size_t, size_t> & size_ima2,
  double * errorMax,
  double * NFA,
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
  std::pair<double,double> ACRansacOut = ACRANSAC(kernel, *pvec_inliers,
    4096, pE, precision, true);
  *errorMax = ACRansacOut.first;
  *NFA = ACRansacOut.second;

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
static bool estimate_Rt_fromE(const Mat3 & K1, const Mat3 & K2,
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
        &x2_ = x2.col(vec_inliers[k]);
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
  if (*iter == 0)
  {
    std::cerr << std::endl << "/!\\There is no right solution,"
      << " probably intermediate results are not correct or no points"
      << " in front of both cameras" << std::endl;
    return false;
  }
  const size_t index = std::distance(f.begin(), iter);
  (*R) = Rs[index];
  (*t) = ts[index];

  return true;
}

} // namespace openMVG

#endif // OPENMVG_SAMPLES_ROBUST_ESSENTIAL_HPP
