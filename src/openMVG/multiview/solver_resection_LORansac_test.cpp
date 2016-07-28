#include <openMVG/robust_estimation/robust_estimator_LORansac.hpp>
#include <openMVG/robust_estimation/robust_estimator_LORansacKernelAdaptor.hpp>
#include <openMVG/robust_estimation/score_evaluator.hpp>
#include <openMVG/robust_estimation/rand_sampling.hpp>
#include <openMVG/multiview/projection.hpp>
#include <openMVG/multiview/solver_resection_kernel.hpp>
#include <openMVG/multiview/solver_resection_p3p.hpp>
#include <openMVG/multiview/conditioning.hpp>
#include <openMVG/cameras/cameras.hpp>
#include <openMVG/sfm/sfm.hpp>
#include <openMVG/geometry/pose3.hpp>
#include <openMVG/numeric/numeric.h>

#include "testing/testing.h"

#include <vector>
#include <random>
#include <algorithm>

using namespace openMVG;

struct ResectionSquaredResidualError
{
  // Compute the residual of the projection distance(pt2D, Project(P,pt3D))
  // Return the squared error

  static double Error(const Mat34 & P, const Vec2 & pt2D, const Vec3 & pt3D)
  {
    const Vec2 x = Project(P, pt3D);
    return (x - pt2D).squaredNorm();
  }
};

bool refinePoseAsItShouldbe(const Mat & pt3D,
                            const Mat & pt2D,
                            const std::vector<std::size_t> & vec_inliers,
                            cameras::IntrinsicBase * intrinsics,
                            geometry::Pose3 & pose,
                            bool b_refine_pose,
                            bool b_refine_intrinsic)
{
  using namespace sfm;

  // Setup a tiny SfM scene with the corresponding 2D-3D data
  SfM_Data sfm_data;
  // view
  sfm_data.views.insert(std::make_pair(0, std::make_shared<View>("", 0, 0, 0)));
  // pose
  sfm_data.poses[0] = pose;
  // intrinsic (the shared_ptr does not take the ownership, will not release the input pointer)
  sfm_data.intrinsics[0] = std::shared_ptr<cameras::IntrinsicBase>(intrinsics, [](cameras::IntrinsicBase*)
  {
  });
  // structure data (2D-3D correspondences)
  for(size_t i = 0; i < vec_inliers.size(); ++i)
  {
    const size_t idx = vec_inliers[i];
    Landmark landmark;
    landmark.X = pt3D.col(idx);
    landmark.obs[0] = Observation(pt2D.col(idx), UndefinedIndexT);
    sfm_data.structure[i] = std::move(landmark);
  }

  Bundle_Adjustment_Ceres bundle_adjustment_obj;
  BA_Refine refineOptions = BA_REFINE_NONE;
  if(b_refine_pose)
    refineOptions |= sfm::BA_REFINE_ROTATION | sfm::BA_REFINE_TRANSLATION;
  if(b_refine_intrinsic)
    refineOptions |= sfm::BA_REFINE_INTRINSICS_ALL;
  const bool b_BA_Status = bundle_adjustment_obj.Adjust(sfm_data, refineOptions);
  if(b_BA_Status)
  {
    pose = sfm_data.poses[0];
  }
  return b_BA_Status;
}


TEST(P3P_Ransac, noisyFromImagePoints)
{
  // camera and image parameters
  const std::size_t WIDTH = 1600;
  const std::size_t HEIGHT = 1200;
  const std::size_t FOCAL = 2000;
  const std::size_t PPX = WIDTH / 2;
  const std::size_t PPY = HEIGHT / 2;
  
  // simulation parameters
  // number of trials to run
  const std::size_t NUMTRIALS = 100;
  // create outliers
  const bool withOutliers = true;
  const double OUTLIERSRATIO = .15;
  assert(OUTLIERSRATIO <= 1.0 && OUTLIERSRATIO >= .0);
  // parameters for generating 3D points
  // 3D points are generating selecting a point on the image and then picking a
  // points on its associated projection ray with a distance in [MINDIST, MAXDIST]
  const double MINDIST = 15;
  const double MAXDIST = 35;
  assert(MINDIST <= MAXDIST);
  // number of points to generate
  const std::size_t nbPoints = 50;
  // noise level in pixels
  const double gaussianNoiseLevel = 4.0;
  // tolerance errors for test to pass
  const double maxAngularError = 0.1;
  const double maxBaselineError = 0.01;
  
  std::vector<std::size_t> vec_outliers;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> realDist(0, 1.0);

  // generate a random RT transform to apply to the 3D points
  const Mat3 Rgt = rotationXYZ(M_PI_2 * realDist(gen), M_PI_2 * realDist(gen), M_PI_2 * realDist(gen));
  const Vec3 Tgt = MINDIST * Vec3(realDist(gen), realDist(gen), realDist(gen));
  Mat3 Kgt;
  Kgt << FOCAL, 0, PPX,
          0, FOCAL, PPY,
          0, 0, 1;
  const Vec3 Cgt = -Rgt.transpose() * Tgt;

  // draw some random points in the camera image
  Mat pts2D = Mat(2, nbPoints);
  Mat pts3D = Mat(3, nbPoints);
  Mat pts3DGt = Mat(3, nbPoints);

  for(std::size_t i = 0; i < nbPoints; ++i)
  {
    // for each point add a 3D point lying on the projection ray and having a 
    // distance [MINDIST, MAXDIST]
    pts2D.col(i) = Vec2(WIDTH * realDist(gen), HEIGHT * realDist(gen));
    Vec3 direction;
    direction << ((pts2D.col(i) - Vec2(PPX, PPY)) / FOCAL), 1;
    direction /= direction.norm();
    direction *= (MINDIST + (MAXDIST - MINDIST) * realDist(gen));
    pts3DGt.col(i) = direction;
    // multiply by the inverse of the pose
    pts3D.col(i) = Rgt.transpose() * direction + -Rgt.transpose() * Tgt;
  }

  // add some gaussian noise to the 2d points
  for(std::size_t i = 0; i < nbPoints; ++i)
  {
    const double theta = realDist(gen)*2 * M_PI;
    const double radius = gaussianNoiseLevel * realDist(gen);
    pts2D.col(i) += radius * Vec2(std::cos(theta), std::sin(theta));
  }
  
  if(withOutliers)
  {
    // take a random sample to be used as outliers
    const std::size_t NUMOUTLIERS = std::size_t(OUTLIERSRATIO*nbPoints);
    robust::UniformSample(NUMOUTLIERS, nbPoints, &vec_outliers);
    std::sort(vec_outliers.begin(), vec_outliers.end());
    
    // add outliers
    for(const auto &idx : vec_outliers)
    {
      std::size_t iter = 0;
      Vec2 pt = Vec2(WIDTH * realDist(gen), HEIGHT * realDist(gen));
      while( (pt-pts2D.col(idx)).norm() < 15*gaussianNoiseLevel)
      {
        pt = Vec2(WIDTH * realDist(gen), HEIGHT * realDist(gen));
        ++iter;
        // safeguard against infinite loops
        if(iter > 1000)
        {
          std::cerr << "Unable to generate a random point, iterations excedeed!";
          assert(false);
        }
      }
      pts2D.col(idx) = pt;
    }
  }

  for(std::size_t trial = 0; trial < NUMTRIALS; ++trial)
  {
    std::cout << "\nTrial #" << trial << std::endl;
    typedef openMVG::euclidean_resection::P3PSolver SolverType;
    typedef openMVG::resection::kernel::SixPointResectionSolver SolverLSType;

    typedef openMVG::robust::KernelAdaptorResectionLORansac_K<SolverType,
                                                              ResectionSquaredResidualError,
                                                              openMVG::robust::UnnormalizerResection,
                                                              SolverLSType,
                                                              Mat34> KernelType;

    // this is just to simplify and use image plane coordinates instead of camera
    // (pixel) coordinates
    Mat pts2Dnorm;
    ApplyTransformationToPoints(pts2D, Kgt.inverse(), &pts2Dnorm);
    KernelType kernel(pts2Dnorm, pts3D, Mat3::Identity());

    std::vector<std::size_t> vec_inliers;
    const double threshold = 2*gaussianNoiseLevel;
    const double normalizedThreshold = Square(threshold / FOCAL);
    robust::ScorerEvaluator<KernelType> scorer(normalizedThreshold);
    Mat34 Pest = robust::LO_RANSAC(kernel, scorer, &vec_inliers);
    
    Mat3 Rest;
    Mat3 Kest;
    Vec3 Test;
    KRt_From_P(Pest, &Kest, &Rest, &Test);

    std::cout << "Est: Pest:\n" << Pest
            << "\nRest:\n" << Rest
            << "\nKest:\n" << Kest
            << "\ntest:\n" << Test << std::endl;
    
    const std::size_t numInliersFound = vec_inliers.size();
    const std::size_t numInliersExpected = nbPoints-vec_outliers.size();
    std::cout << "Solution found with " << numInliersFound << " inliers" << std::endl;
    std::cout << "Expected number of inliers " << numInliersExpected << std::endl;

    CHECK_EQUAL(numInliersFound, numInliersExpected);
    
    const double angError = R2D(getRotationMagnitude(Rgt * Rest.transpose()));
    std::cout << "Angular error: " << angError << std::endl
            << "baseline error: " << (Test - Tgt).squaredNorm() << std::endl;

    geometry::Pose3 pose = geometry::poseFromRT(Rest, Test);
    refinePoseAsItShouldbe(pts3D,
                           pts2Dnorm,
                           vec_inliers,
                           new cameras::Pinhole_Intrinsic(WIDTH, HEIGHT, 1, 0, 0),
                           pose,
                           true,
                           false );

    const double angErrorRef = R2D(getRotationMagnitude(Rgt * pose.rotation().transpose()));
    const double baselineErrorRef = (Tgt - pose.translation()).squaredNorm();
    std::cout << "Final angular error #"<<trial<<" : " << angErrorRef << std::endl
            << "Final baseline error #"<<trial<<" : " << baselineErrorRef << std::endl;
    
    EXPECT_NEAR(angErrorRef, 0.0, maxAngularError);
    EXPECT_NEAR(baselineErrorRef, 0.0, maxBaselineError);

    std::cout << "Refined pose:\n"
                << "\nEst: Rest:\n" << pose.rotation()
                << "\nCest:\n" << pose.center()
                << "\nTest:\n" << pose.translation() << std::endl;
    
    if(withOutliers)
    {
      // test if inliers found and outliers GT have a empty intersection
      std::vector<std::size_t> inters(nbPoints);
      std::sort(vec_inliers.begin(), vec_inliers.end());
      auto it = std::set_intersection(vec_inliers.begin(), vec_inliers.end(),
                                      vec_outliers.begin(), vec_outliers.end(),
                                      inters.begin());
      inters.resize(it-inters.begin());
      if(inters.size()>0)
        std::cerr << "******* there are " << inters.size() << " outliers considered as inliers" << std::endl;
      CHECK_EQUAL(inters.size(), 0);
    }
  }
}

/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */
