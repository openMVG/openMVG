
// Copyright (c) 2012, 2013 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <iterator>

#include "openMVG/robust_estimation/robust_estimator_lineKernel_test.hpp"
#include "openMVG/robust_estimation/robust_estimator_ACRansac.hpp"
#include "testing/testing.h"

#include "third_party/vectorGraphics/svgDrawer.hpp"
using namespace svg;

using namespace openMVG;
using namespace openMVG::robust;
using namespace std;

/// ACRansac Kernel for line estimation
template <typename SolverArg,
  typename ErrorArg,
  typename ModelArg >
class ACRANSACOneViewKernel
{
public:
  typedef SolverArg Solver;
  typedef ModelArg  Model;

  ACRANSACOneViewKernel(const Mat &x1, int w1, int h1)
    : x1_(x1.rows(), x1.cols()), N1_(3,3), logalpha0_(0.0)
  {
    assert(2 == x1_.rows());

    x1_ = x1;
    N1_ = Mat3::Identity();

    // Model error as point to line error
    // Ratio of containing diagonal image rectangle over image area
    double D = sqrt(w1 *w1*1.0 + h1 * h1); // diameter
    double A = w1 * h1; // area
    logalpha0_ = log10(2.0*D/A /1.0);
  }

  enum { MINIMUM_SAMPLES = Solver::MINIMUM_SAMPLES };
  enum { MAX_MODELS = Solver::MAX_MODELS };

  void Fit(const std::vector<size_t> &samples, std::vector<Model> *models) const {
    Mat sampled_xs = ExtractColumns(x1_, samples);
    Solver::Solve(sampled_xs, models);
  }

  double Error(size_t sample, const Model &model) const {
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
  std::vector<size_t> vec_inliers;
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
  std::vector<size_t> vec_inliers;
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
  std::vector<size_t> vec_inliers;
  Vec2 line;
  ACRANSAC(lineKernel, vec_inliers, 300, &line);

  EXPECT_NEAR(2.0, line[1], 1e-9);
  EXPECT_NEAR(1.0, line[0], 1e-9);
  CHECK_EQUAL(5, vec_inliers.size());
}


// Test if the robust estimator do not return inlier if too few point
// was given for an estimation.
TEST(RansacLineFitter, TooFewPoints) {

  Mat2X xy(2, 1);
  // y = 2x + 1
  xy << 1,
    3;
  ACRANSACOneViewKernel<LineSolver, pointToLineError, Vec2> lineKernel(xy, 12, 12);
  std::vector<size_t> vec_inliers;
  ACRANSAC(lineKernel, vec_inliers);

  CHECK_EQUAL(0, vec_inliers.size());
}

// From a GT model :
//  Compute a list of point that fit the model.
//  Add white noise to given amount of points in this list.
//  Check that the number of inliers and the model are correct.
TEST(RansacLineFitter, RealisticCase) {

  const int NbPoints = 100;
  const int inlierPourcentAmount = 30;
  Mat2X xy(2, NbPoints);

  Vec2 GTModel;
  GTModel <<  -2.0, 6.3;

  //-- Build the point list according the given model
  for(int i = 0; i < NbPoints; ++i) {
    xy.col(i) << i, (double)i*GTModel[1] + GTModel[0];
  }

  //-- Add some noise (for the asked percentage amount)
  int nbPtToNoise = (int) NbPoints*inlierPourcentAmount/100.0;
  vector<size_t> vec_samples; // Fit with unique random index
  UniformSample(nbPtToNoise, NbPoints, &vec_samples);
  for(size_t i = 0; i <vec_samples.size(); ++i)
  {
    const size_t randomIndex = vec_samples[i];
    //Additive random noise
    xy.col(randomIndex) += Vec2::Random()/10.;
  }

  // The base estimator
  ACRANSACOneViewKernel<LineSolver, pointToLineError, Vec2> lineKernel(xy, 12, 12);

  // Check the best model that fit the most of the data
  //  in a robust framework (ACRANSAC).
  std::vector<size_t> vec_inliers;
  Vec2 line;
  ACRANSAC(lineKernel, vec_inliers, 300, &line);

  CHECK_EQUAL(NbPoints-nbPtToNoise, vec_inliers.size());
  EXPECT_NEAR(GTModel(0), line[0], 1e-9);
  EXPECT_NEAR(GTModel(1), line[1], 1e-9);
}

// Generate a random value between [0;1]
static inline double randValue()  {
  return rand()/double(RAND_MAX);
}

double gaussGenerator(){
  double f_x1, f_x2, f_w;
  static int i_set = 0;
  static float f_noise;

  if(i_set==0){// We do not have an extra deviate handy so,
    do {
      f_x1 = 2.0 * randValue() - 1.0;// pick 2 uniform numbers in the square
      f_x2 = 2.0 * randValue() - 1.0;// extending from -1 to 1 in each direction
      f_w = f_x1 * f_x1 + f_x2 * f_x2;// see if they are in the unit circle
    } while ( (f_w >= 1.0) || (f_w == 0.0));// and if they are not, try again

    f_w = sqrt( (-2.0 * log( f_w ) ) / f_w );

    // Now make the Box-Muller transformation to get 2 normal deviates
    f_noise=f_x2*f_w;// save one for the next  time
    i_set=1;// set flag
    return (f_x1*f_w);// return one
  }else{// We have an extra deviate handy,
    i_set=0;// so unset the flag
    return (f_noise);// and return the value
  }
}

double generateGaussianErrorAroundAValue(double v, double d_deviation)
{
  double d_value, d_noise;
  d_value = (double) v;
  d_noise = gaussGenerator();
  d_noise *= d_deviation;
  d_value += d_noise;
  return(d_value);
}

// Generate nbPoints along a line and add gaussian noise.
// Move some point in the dataset to create outlier contamined data
void generateLine(Mat & points, size_t nbPoints, int W, int H, float noise, float outlierPourcent)
{
  points = Mat(2, nbPoints);

  Vec2 lineEq(50,0.3);

  for (size_t i = 0; i < nbPoints; ++i)
  {
    float x = rand()%W;
    float y =  lineEq[1] * x + lineEq[0];
    y = generateGaussianErrorAroundAValue(y, noise);
    points.col(i) = Vec2(x, y);
  }

  // generate outlier
  size_t count = outlierPourcent * nbPoints;
  std::vector<size_t> vec_indexes(count,0);
  random_sample(count, nbPoints, &vec_indexes);
  for (size_t i = 0; i < count; ++i)
  {
    size_t pos = vec_indexes[i];
    float noise2 = generateGaussianErrorAroundAValue(0, 0.2);
    float noise3 = generateGaussianErrorAroundAValue(0, 0.2);
    points.col(pos) = Vec2(rand()%W +noise2, rand()%H - noise3);
  }
}

// Structure used to avoid repetition in a given series
struct IndMatchd
{
  IndMatchd(double i = 0, double j = 0): _i(i), _j(j)
  {}

  friend bool operator==(const IndMatchd& m1, const IndMatchd& m2)
  {    return (m1._i == m2._i && m1._j == m2._j);  }

  // Lexicographical ordering of matches. Used to remove duplicates.
  friend bool operator<(const IndMatchd& m1, const IndMatchd& m2)
  {
    if(m1._i < m2._i) return true;
    if(m1._i > m2._i) return false;

    if(m1._j < m2._j) return true;
    else
      return false;

  }

  double _i, _j;
};

// Test ACRANSAC adaptability to noise
// Set a line with a increasing gaussian noise
// See if the AContrario RANSAC is able to label the good point as inlier
//  by having it's estimated confidence threshold growing.
TEST(RansacLineFitter, ACRANSACSimu) {

  int S = 100;
  int W = S, H = S;

  std::vector<double> vec_gaussianValue;
  for (int i=0; i < 10; ++i)  {
    vec_gaussianValue.push_back(i/10. * 5.);
  }

  for (std::vector<double>::const_iterator iter = vec_gaussianValue.begin();
    iter != vec_gaussianValue.end(); ++ iter)
  {
    double gaussianNoiseLevel = *iter;

    size_t nbPoints = 2.0 * S * sqrt(2.0);
    float noise = gaussianNoiseLevel;

    float outlierPourcent = .3;
    Mat points;
    generateLine(points, nbPoints, W, H, noise, outlierPourcent);

    // Remove point that have the same coords
    {
      std::vector<IndMatchd> vec_match;
      for (size_t i = 0; i < nbPoints; ++i) {
        vec_match.push_back(IndMatchd(points.col(i)[0], points.col(i)[1]));
      }

      std::sort(vec_match.begin(), vec_match.end());
      std::vector<IndMatchd>::iterator end = std::unique(vec_match.begin(), vec_match.end());
      if(end != vec_match.end()) {
        std::cout << "Remove " << std::distance(end, vec_match.end())
          << "/" << vec_match.size() << " duplicate matches, "
          << " keeping " << std::distance(vec_match.begin(), end) <<std::endl;
        vec_match.erase(end, vec_match.end());
      }

      points = Mat(2, vec_match.size());
      for (size_t i = 0; i  < vec_match.size(); ++i)  {
        points.col(i) = Vec2(vec_match[i]._i, vec_match[i]._j);
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
    std::vector<size_t> vec_inliers;
    std::pair<double,double> ret = ACRANSAC(lineKernel, vec_inliers, 1000, &line);
    double errorMax = ret.first;

    EXPECT_TRUE( abs(gaussianNoiseLevel - sqrt(errorMax)) < 2.0);

    ostringstream os;
    os << gaussianNoiseLevel << "_line_" << sqrt(errorMax) << ".png";

    //Svg drawing
    {
      svgDrawer svgTest(W,H);
      for (size_t i = 0; i < nbPoints; ++i) {
        string sCol = "red";
        float x = points.col(i)[0];
        float y = points.col(i)[1];
        if (find(vec_inliers.begin(), vec_inliers.end(), i) != vec_inliers.end())
        {
          sCol = "green";
        }
        svgTest.drawCircle(x,y, 1, svgStyle().fill(sCol).noStroke());
      }
      //draw the found line
      float xa = 0, xb = W;
      float ya = line[1] * xa + line[0];
      float yb = line[1] * xb + line[0];
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
