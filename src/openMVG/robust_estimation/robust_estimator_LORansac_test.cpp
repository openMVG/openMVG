#include "openMVG/robust_estimation/robust_estimator_lineKernel_test.hpp"
#include "openMVG/robust_estimation/robust_estimator_LORansac.hpp"
#include "openMVG/robust_estimation/score_evaluator.hpp"

#include "openMVG/numeric/numeric.h"
#include "third_party/vectorGraphics/svgDrawer.hpp"
#include "testing/testing.h"

#include <iostream>
#include <random>
#include <fstream>
#include <vector>
#include <string>

using namespace openMVG;
using namespace openMVG::robust;

/**
 * @brief Generate a svg file with the ground truth line, the estimated one, the
 * estimated inliers and outliers.
 * 
 * @param[in] outfile The name of the svg file to generate.
 * @param[in] W The width of the image to generate.
 * @param[in] H The height of the image to generate.
 * @param[in] lineGT The ground truth line.
 * @param[in] lineEst The estimated line.
 * @param[in] points The points from which the lines are generated.
 * @param[in] vec_inliers The inliers that fit the estimated line.
 */
void drawTest(const std::string &outfile,
              int imageWidth, 
              int imageHeight,
              const Vec2 &lineGT,
              const Vec2 &lineEst,
              const Mat &points,
              const std::vector<std::size_t> &vec_inliers)
{
  const std::size_t nbPoints = points.cols();
  svg::svgDrawer svgTest(imageWidth, imageHeight);
  for(std::size_t i = 0; i < nbPoints; ++i)
  {
    std::string sCol = "red";
    float x = points.col(i)[0];
    float y = points.col(i)[1];
    if(std::find(vec_inliers.begin(), vec_inliers.end(), i) != vec_inliers.end())
    {
      sCol = "green";
    }
    svgTest.drawCircle(x, y, 1, svg::svgStyle().fill(sCol).noStroke());
  }
  //draw the found line
  float xa = 0, xb = imageWidth;
  float ya = lineEst[1] * xa + lineEst[0];
  float yb = lineEst[1] * xb + lineEst[0];
  svgTest.drawLine(xa, ya, xb, yb, svg::svgStyle().stroke("blue", 0.5));
  //draw the GT line
  ya = lineGT[1] * xa + lineGT[0];
  yb = lineGT[1] * xb + lineGT[0];
  svgTest.drawLine(xa, ya, xb, yb, svg::svgStyle().stroke("black", 0.5));

  //  ostringstream osSvg;
  //  osSvg << gaussianNoiseLevel << "_line_" << sqrt(errorMax) << ".svg";
  std::ofstream svgFile(outfile);
  svgFile << svgTest.closeSvgFile().str();
  svgFile.close();
}

struct LineKernelLoRansac : public LineKernel
{
  typedef Vec2 Model; // line parametrization: a, b;

  enum
  {
    MINIMUM_SAMPLES = 2,
    MINIMUM_LSSAMPLES = 2
  };

  LineKernelLoRansac(const Mat2X &xs) : LineKernel(xs)
  {
  }

  void FitLS(const std::vector<size_t> &samples, std::vector<Vec2> *lines, const std::vector<double> *weights = nullptr) const
  {

    assert(samples.size() >= (unsigned int) MINIMUM_SAMPLES);
    // Standard least squares solution.
    const Mat2X sampled_xs = ExtractColumns(xs_, samples);
    if(weights)
      LineSolver::SolveWeightedLS(sampled_xs, lines, *weights);
    else
      LineSolver::Solve(sampled_xs, lines);
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
      vec_weights[sample] = Error(idx, model);
      // avoid division by zero
      vec_weights[sample] = 1 / std::pow(std::max(eps, vec_weights[sample]), 2);
    }
  }
};

TEST(LoRansacLineFitter, IdealCaseLoRansac) 
{

  const int NbPoints = 3000;
  const int outlierRatio = 30; //works with 40
  Mat2X xy(2, NbPoints);

  Vec2 GTModel; // y = 2x + 6.3
  GTModel <<  -2.0, 6.3;

  //-- Build the point list according the given model
  for(std::size_t i = 0; i < NbPoints; ++i)  
  {
    xy.col(i) << i, (double)i*GTModel[1] + GTModel[0];
  }

  //-- Add some noise (for the asked percentage amount)
  int nbPtToNoise = (int) NbPoints*outlierRatio/100.0;
  vector<size_t> vec_samples; // Fit with unique random index
  UniformSample(nbPtToNoise, NbPoints, &vec_samples);
  for(size_t i = 0; i <vec_samples.size(); ++i)
  {
    const size_t randomIndex = vec_samples[i];
    //Additive random noise
    xy.col(randomIndex) << xy.col(randomIndex)(0)+rand()%2-3,
                           xy.col(randomIndex)(1)+rand()%8-6;
  }
  
  LineKernelLoRansac kernel(xy);
  std::vector<size_t> vec_inliers;
  Vec2 model = LO_RANSAC(kernel, ScorerEvaluator<LineKernel>(0.3), &vec_inliers);
  std::cout << "#inliers found : " << vec_inliers.size() 
          << " expected: " << NbPoints-nbPtToNoise << std::endl;
  std::cout << "model[0] found : " << model[0] 
          << " expected: " << GTModel[0] << std::endl;
  std::cout << "model[1] found : " << model[1] 
          << " expected: " << GTModel[1] << std::endl;
  
  CHECK_EQUAL(NbPoints-nbPtToNoise, vec_inliers.size());
  EXPECT_NEAR(GTModel[0], model[0], 1e-9);
  EXPECT_NEAR(GTModel[1], model[1], 1e-9);
}

TEST(LoRansacLineFitter, RealCaseLoRansac)
{

  const int NbPoints = 3000;
  const int outlierRatio = 30;
  const double gaussianNoiseLevel = 0.1;
  const std::size_t numTrials = 1;
  
  Mat2X xy(2, NbPoints);

  Vec2 GTModel; // y = 2x + 1
  GTModel << -2, .3;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(0, gaussianNoiseLevel);
  std::uniform_real_distribution<double> realDist(0, 1.0);
  
  for(std::size_t trial = 0; trial < numTrials; ++trial)
  {

    //-- Build the point list according the given model adding some noise
    for(std::size_t i = 0; i < NbPoints; ++i)
    {
      xy.col(i) << i, (double) i * GTModel[1] + GTModel[0];
      const double theta = realDist(gen)*2 * M_PI;
      const double radius = d(gen);
      xy.col(i) += radius * Vec2(std::cos(theta), std::sin(theta));
    }
    const int W = std::abs(xy(0, 0) - xy(0, NbPoints - 1));
    const int H = (int) std::fabs(xy(1, 0) - xy(1, NbPoints - 1));
    ;

    //-- Add some outliers (for the asked percentage amount)
    int nbPtToNoise = (int) NbPoints * outlierRatio / 100.0;
    vector<std::size_t> vec_inliersGT(NbPoints);
    std::iota(vec_inliersGT.begin(), vec_inliersGT.end(), 0);

    vector<std::size_t> vec_outliers; // Fit with unique random index
    UniformSample(nbPtToNoise, NbPoints, &vec_outliers);

    for(std::size_t i = 0; i < vec_outliers.size(); ++i)
    {
      const std::size_t randomIndex = vec_outliers[i];

      vec_inliersGT.erase(std::remove(vec_inliersGT.begin(), vec_inliersGT.end(), randomIndex), vec_inliersGT.end());

      Vec2 pt;
      double distance = 0;
      // try to generate a point that is well far from the line
      while(distance < 15 * gaussianNoiseLevel)
      {
        pt(0) = realDist(gen) * W;
        pt(1) = realDist(gen) * H;
        distance = pointToLineError::Error(GTModel, pt);
      }

      xy.col(randomIndex) = pt;
    }
    // just sanity check
    CHECK_EQUAL(NbPoints - nbPtToNoise, vec_inliersGT.size());

    std::map<std::size_t, std::size_t> histogram;
    for(std::size_t i = 0; i < NbPoints; ++i)
      histogram[i] = 0;

    const double threshold = 3 * gaussianNoiseLevel;
    const std::string base = "testRansac_line_t" + std::to_string(threshold) + "_n" + std::to_string(gaussianNoiseLevel);
    LineKernelLoRansac kernel(xy);
    std::vector<size_t> vec_inliers;

    Vec2 model = LO_RANSAC(kernel, ScorerEvaluator<LineKernel>(threshold), &vec_inliers);
    std::cout << "#inliers found : " << vec_inliers.size()
            << " expected: " << NbPoints - nbPtToNoise << std::endl;
    std::cout << "model[0] found : " << model[0]
            << " expected: " << GTModel[0] << std::endl;
    std::cout << "model[1] found : " << model[1]
            << " expected: " << GTModel[1] << std::endl;

    drawTest(base + "_LORANSACtrial" + std::to_string(0) + ".svg",
             W, H,
             GTModel,
             model,
             xy,
             vec_inliers);


    CHECK_EQUAL(NbPoints - nbPtToNoise, vec_inliers.size());
  }
}

/* ************************************************************************* */
int main()
{
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
/* ************************************************************************* */