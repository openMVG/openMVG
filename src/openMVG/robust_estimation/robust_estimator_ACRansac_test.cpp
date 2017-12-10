// This file is part of OpenMVG, an Open Multiple View Geometry C++ library.

// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "openMVG/robust_estimation/robust_estimator_lineKernel_test.hpp"
#include "openMVG/numeric/extract_columns.hpp"

#include "testing/testing.h"
#include "third_party/vectorGraphics/svgDrawer.hpp"

#include <iterator>
#include <random>


using namespace openMVG;
using namespace openMVG::robust;
using namespace std;
using namespace svg;

/// ACRansac Kernel for line estimation
template <typename SolverArg,
  typename ErrorArg,
  typename ModelArg >
class ACRANSACOneViewKernel
{
public:
  using Solver = SolverArg;
  using Model = ModelArg;

  ACRANSACOneViewKernel(const Mat &x1, int w1, int h1)
    : x1_(x1), N1_(Mat3::Identity()), logalpha0_(0.0)
  {
    assert(2 == x1_.rows());

    // Model error as point to line error
    // Ratio of containing diagonal image rectangle over image area
    const double D = sqrt(w1 *w1*1.0 + h1 * h1); // diameter
    const double A = w1 * h1; // area
    logalpha0_ = log10(2.0*D/A /1.0);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<uint32_t> &samples, std::vector<Model> *models) const {
    const Mat sampled_xs = ExtractColumns(x1_, samples);
    Solver::Solve(sampled_xs, models);
  }

  double Error(uint32_t sample, const Model &model) const {
    return ErrorArg::Error(model, x1_.col(sample));
  }

  void Errors(const Model &model, std::vector<double> & vec_errors) const {
    for (size_t sample = 0; sample < x1_.cols(); ++sample)
      vec_errors[sample] = ErrorArg::Error(model, x1_.col(sample));
  }

  size_t NumSamples() const {
    return x1_.cols();
  }

  void Unnormalize(Model * model) const {
    // Model is left unchanged
  }

  double logalpha0() const {return logalpha0_;}

  double multError() const {return 0.5;}

  Mat3 normalizer1() const {return Mat3::Identity();}
  Mat3 normalizer2() const {return Mat3::Identity();}
  double unormalizeError(double val) const {return sqrt(val);}

private:
  Mat x1_;
  Mat3 N1_;
  double logalpha0_;
};

// Test ACRANSAC with the AC-adapted Line kernel in a noise/outlier free dataset
TEST(RansacLineFitter, OutlierFree) {

  Mat2X xy(2, 5);
  // y = 2x + 1
  xy << 1, 2, 3, 4,  5,
        3, 5, 7, 9, 11;

  // The base estimator
  ACRANSACOneViewKernel<LineSolver, pointToLineError, Vec2> lineKernel(xy, 12, 12);

  // Check the best model that fit the most of the data
  //  in a robust framework (ACRANSAC).
  std::vector<uint32_t> vec_inliers;
  Vec2 line;
  ACRANSAC(lineKernel, vec_inliers, 300, &line);

  EXPECT_NEAR(2.0, line[1], 1e-9);
  EXPECT_NEAR(1.0, line[0], 1e-9);
  CHECK_EQUAL(5, vec_inliers.size());
}

// Simple test without getting back the model
TEST(RansacLineFitter, OutlierFree_DoNotGetBackModel) {

  Mat2X xy(2, 5);
  // y = 2x + 1
  xy << 1, 2, 3, 4,  5,
        3, 5, 7, 9, 11;

  ACRANSACOneViewKernel<LineSolver, pointToLineError, Vec2> lineKernel(xy, 12, 12);
  std::vector<uint32_t> vec_inliers;
  ACRANSAC(lineKernel, vec_inliers);

  CHECK_EQUAL(5, vec_inliers.size());
}


TEST(RansacLineFitter, OneOutlier) {

  Mat2X xy(2, 6);
  // y = 2x + 1 with an outlier
  xy << 1, 2, 3, 4,  5, 100, // (100,-123) is the outlier
        3, 5, 7, 9, 11, -123;

  // The base estimator
  ACRANSACOneViewKernel<LineSolver, pointToLineError, Vec2> lineKernel(xy, 12, 12);

  // Check the best model that fit the most of the data
  //  in a robust framework (ACRANSAC).
  std::vector<uint32_t> vec_inliers;
  Vec2 line;
  ACRANSAC(lineKernel, vec_inliers, 300, &line);

  EXPECT_NEAR(2.0, line[1], 1e-9);
  EXPECT_NEAR(1.0, line[0], 1e-9);
  CHECK_EQUAL(5, vec_inliers.size());
}


// Test if the robust estimator do not return inlier if too few point
// was given for an estimation.
TEST(RansacLineFitter, TooFewPoints) {

  Vec2 xy;
  // y = 2x + 1
  xy << 1, 2;
  ACRANSACOneViewKernel<LineSolver, pointToLineError, Vec2> lineKernel(xy, 12, 12);
  std::vector<uint32_t> vec_inliers;
  ACRANSAC(lineKernel, vec_inliers);

  CHECK_EQUAL(0, vec_inliers.size());
}

// From a GT model :
//  Compute a list of point that fit the model.
//  Add white noise to given amount of points in this list.
//  Check that the number of inliers and the model are correct.
TEST(RansacLineFitter, RealisticCase) {

  constexpr int NbPoints = 100;
  constexpr double inlierRatio = 30.0 / 100.0;
  Mat2X xy(2, NbPoints);

  Vec2 GTModel; // y = 6.3 x + (-2.0)
  GTModel <<  -2.0, 6.3;

  //-- Build the point list according the given model
  for (int i = 0; i < NbPoints; ++i) {
    xy.col(i) << i, static_cast<double>(i)*GTModel[1] + GTModel[0];
  }

  // Setup a normal distribution in order to make outlier not aligned
  std::mt19937 random_generator(std::mt19937::default_seed);
  std::normal_distribution<> d(0, 5); // More or less 5 units

  //-- Simulate outliers (for the asked percentage amount of the datum)
  constexpr auto nbPtToNoise = static_cast<uint32_t>(NbPoints*inlierRatio);
  vector<uint32_t> vec_samples; // Fit with unique random index
  UniformSample(nbPtToNoise, NbPoints, random_generator, &vec_samples);
  for (size_t i = 0; i <vec_samples.size(); ++i)
  {
    const size_t randomIndex = vec_samples[i];
    // Start from a outlier point (0,0)
    // and move it in a given small range (since it must remains in an outlier area)
    xy.col(randomIndex)<< d(random_generator), d(random_generator);
  }

  // The base estimator
  ACRANSACOneViewKernel<LineSolver, pointToLineError, Vec2> lineKernel(xy, 12, 12);

  // Check the best model that fit the most of the data
  //  in a robust framework (ACRANSAC).
  std::vector<uint32_t> vec_inliers;
  Vec2 line;
  ACRANSAC(lineKernel, vec_inliers, 300, &line);

  CHECK_EQUAL(NbPoints-nbPtToNoise, vec_inliers.size());
  EXPECT_NEAR(GTModel(0), line[0], 1e-9);
  EXPECT_NEAR(GTModel(1), line[1], 1e-9);
}

// Generate nbPoints along a line and add gaussian noise.
// Move some point in the dataset to create outlier contamined data
void generateLine(Mat & points, size_t nbPoints, int W, int H, float noise, float outlierRatio)
{
  points = Mat(2, nbPoints);

  Vec2 lineEq(50, 0.3);

  // Setup a normal distribution of mean 0 and amplitude equal to noise
  std::mt19937 random_generator(std::mt19937::default_seed);
  std::normal_distribution<double> d(0, noise);

  // Setup uniform distribution
  std::uniform_int_distribution<int> dW(0, W), dH(0, H);

  for (size_t i = 0; i < nbPoints; ++i)
  {
    auto x = static_cast<double>(dW(random_generator));
    const float y =  d(random_generator) + (lineEq[1] * x + lineEq[0]) + d(random_generator);
    points.col(i) << x, y;
  }

  // generate outlier
  std::normal_distribution<double> d_outlier(0, 0.2);
  const auto count = static_cast<uint32_t>(outlierRatio * nbPoints);
  std::vector<uint32_t> vec_indexes(count,0);
  UniformSample(count, nbPoints, random_generator, &vec_indexes);
  for (const auto & pos : vec_indexes)
  {
    points.col(pos) << dW(random_generator) + d_outlier(random_generator),
                       dH(random_generator) - d_outlier(random_generator);
  }
}

// Structure used to avoid repetition in a given series
struct IndMatchd
{
  IndMatchd(double i = 0, double j = 0): i_(i), j_(j)
  {}

  friend bool operator==(const IndMatchd& m1, const IndMatchd& m2)
  {    return (m1.i_ == m2.i_ && m1.j_ == m2.j_);  }

  // Lexicographical ordering of matches. Used to remove duplicates.
  friend bool operator<(const IndMatchd& m1, const IndMatchd& m2)
  {
    if (m1.i_ < m2.i_) return true;
    if (m1.i_ > m2.i_) return false;
    return (m1.j_ < m2.j_);
  }

  double i_, j_;
};

// Test ACRANSAC adaptability to noise
// Set a line with a increasing gaussian noise
// See if the AContrario RANSAC is able to label the good point as inlier
//  by having its estimated confidence threshold growing.
TEST(RansacLineFitter, ACRANSACSimu) {

  const int S = 100;
  const int W = S, H = S;

  std::vector<double> vec_gaussianValue;
  for (int i=0; i < 10; ++i)  {
    vec_gaussianValue.push_back(i/10. * 5. + std::numeric_limits<double>::epsilon());
  }

  for (std::vector<double>::const_iterator iter = vec_gaussianValue.begin();
    iter != vec_gaussianValue.end(); ++ iter)
  {
    const double gaussianNoiseLevel = *iter;

    size_t nbPoints = 2.0 * S * sqrt(2.0);
    const float noise = gaussianNoiseLevel;

    const float outlierRatio = .3f;
    Mat points;
    generateLine(points, nbPoints, W, H, noise, outlierRatio);

    // Remove point that have the same coords
    {
      std::vector<IndMatchd> vec_match;
      for (size_t i = 0; i < nbPoints; ++i) {
        vec_match.push_back(IndMatchd(points.col(i)[0], points.col(i)[1]));
      }

      std::sort(vec_match.begin(), vec_match.end());
      std::vector<IndMatchd>::iterator end = std::unique(vec_match.begin(), vec_match.end());
      if (end != vec_match.end()) {
        std::cout << "Remove " << std::distance(end, vec_match.end())
          << "/" << vec_match.size() << " duplicate matches, "
          << " keeping " << std::distance(vec_match.begin(), end) <<std::endl;
        vec_match.erase(end, vec_match.end());
      }

      points = Mat(2, vec_match.size());
      for (size_t i = 0; i  < vec_match.size(); ++i)  {
        points.col(i) = Vec2(vec_match[i].i_, vec_match[i].j_);
      }
      nbPoints = vec_match.size();
    }

    // draw image
    svgDrawer svgTest(W,H);
    for (size_t i = 0; i < nbPoints; ++i) {
      float x = points.col(i)[0], y = points.col(i)[1];
      svgTest.drawCircle(x,y, 1, svgStyle().fill("red").noStroke());
    }

    ostringstream osSvg;
    osSvg << gaussianNoiseLevel << "_line_.svg";
    ofstream svgFile( osSvg.str().c_str());
    svgFile << svgTest.closeSvgFile().str();

    // robust line estimation
    Vec2 line;

    // The base estimator
    ACRANSACOneViewKernel<LineSolver, pointToLineError, Vec2> lineKernel(points, W, H);

    // Check the best model that fit the most of the data
    //  in a robust framework (ACRANSAC).
    std::vector<uint32_t> vec_inliers;
    const std::pair<double,double> ret = ACRANSAC(lineKernel, vec_inliers, 1000, &line);
    const double errorMax = ret.first;
    EXPECT_TRUE(sqrt(errorMax) < noise*2+.5);

    ostringstream os;
    os << gaussianNoiseLevel << "_line_" << sqrt(errorMax) << ".png";

    // Svg drawing
    {
      svgDrawer svgTest(W,H);
      for (size_t i = 0; i < nbPoints; ++i) {
        string sCol = "red";
        const float x = points.col(i)[0];
        const float y = points.col(i)[1];
        if (find(vec_inliers.begin(), vec_inliers.end(), i) != vec_inliers.end())
        {
          sCol = "green";
        }
        svgTest.drawCircle(x,y, 1, svgStyle().fill(sCol).noStroke());
      }
      //draw the found line
      const float xa = 0, xb = W;
      const float ya = line[1] * xa + line[0];
      const float yb = line[1] * xb + line[0];
      svgTest.drawLine(xa, ya, xb, yb, svgStyle().stroke("blue", 0.5));

      ostringstream osSvg;
      osSvg << gaussianNoiseLevel << "_line_" << sqrt(errorMax) << ".svg";
      ofstream svgFile( osSvg.str().c_str());
      svgFile << svgTest.closeSvgFile().str();
    }
  }
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
