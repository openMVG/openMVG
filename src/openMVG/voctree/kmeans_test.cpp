/* 
 * File:   kmean_test.cpp
 * Author: sgaspari
 *
 * Created on February 11, 2016, 11:54 AM
 */

#include <openMVG/voctree/simple_kmeans.hpp>

#include <testing/testing.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <random>

TEST(kmeans, kmeanInitializer) 
{
  using namespace openMVG;
  
  std::cout << "Testing kmeanspp Initializer..." << std::endl;

  const int DIMENSION = 128;
  const int FEATURENUMBER = 1000;

  const uint32_t K = 10;

  typedef Eigen::Matrix<float, 1, DIMENSION> FeatureFloat;
  typedef std::vector<FeatureFloat, Eigen::aligned_allocator<FeatureFloat> > FeatureFloatVector;
  typedef std::vector<FeatureFloat* > FeaturePointerVector;

  FeatureFloatVector features;
  FeaturePointerVector featPtr;
  FeatureFloatVector centers;
  features.reserve(FEATURENUMBER);
  featPtr.reserve(features.size());
  
  for(std::size_t i = 0; i < FEATURENUMBER; ++i)
  {
    features.push_back(FeatureFloat::Random(1, DIMENSION));
    featPtr.push_back(const_cast<FeatureFloat*> (&features.back()));
  }

  voctree::InitKmeanspp initializer;

  initializer(featPtr, K, centers, voctree::L2<FeatureFloat, FeatureFloat>());

  // it's difficult to check the result as it is random, just check there are no weird things
  EXPECT_TRUE(voctree::checkVectorElements(centers, "initializer1"));

  // now try to generate k cluster well far away and comapare
  features.clear();
  featPtr.clear();
  features.reserve(FEATURENUMBER * K);
  featPtr.reserve(features.size());
  for(int i = 0; i < K; ++i)
  {
    // at each i iteration translate the cluster by 5*i
    for(int j = 0; j < FEATURENUMBER; ++j)
    {
      features.push_back(FeatureFloat::Random(1, DIMENSION) + Eigen::MatrixXf::Constant(1, DIMENSION, 5 * i));
      featPtr.push_back(const_cast<FeatureFloat*> (&features.back()));
    }
  }

  initializer(featPtr, K, centers, voctree::L2<FeatureFloat,FeatureFloat>());

  // it's difficult to check the result as it is random, just check there are no weird things
  EXPECT_TRUE(voctree::checkVectorElements(centers, "initializer2"));
  
}


TEST(kmeans, kmeanInitializerVarying)
{
  using namespace openMVG;
  
  std::cout << "Testing kmeanspp Initializer with variable k and DIM..." << std::endl;

  const int FEATURENUMBER = 1000;
  const std::size_t numTrial = 10;
  using namespace std;

  // generate random values for K and DIMENSION
  std::default_random_engine generator;
  std::uniform_int_distribution<std::size_t> dimGen(3, 128);
  std::uniform_int_distribution<std::size_t> kGen(6, 300);

  for(size_t trial = 0; trial < numTrial; ++trial)
  {
    const std::size_t DIMENSION = dimGen(generator);
    const std::size_t K = kGen(generator);
    const std::size_t STEP = 5 * K;
    std::cout << "\tTrial " << trial + 1 << "/" << numTrial << " with K = " << K << " and DIMENSION = " << DIMENSION << std::endl;

    typedef Eigen::RowVectorXf FeatureFloat;
    typedef std::vector<FeatureFloat, Eigen::aligned_allocator<FeatureFloat> > FeatureFloatVector;
    typedef std::vector<FeatureFloat* > FeaturePointerVector;

    FeatureFloatVector features;
    FeaturePointerVector featPtr;
    FeatureFloatVector centers;

    voctree::InitKmeanspp initializer;

    features.reserve(FEATURENUMBER * K);
    featPtr.reserve(features.size());
    featPtr.reserve(features.size());
    for(std::size_t i = 0; i < K; ++i)
    {
      // at each i iteration translate the cluster by 5*i
      for(std::size_t j = 0; j < FEATURENUMBER; ++j)
      {
        features.push_back((FeatureFloat::Random(DIMENSION) + FeatureFloat::Constant(DIMENSION, STEP * i) - FeatureFloat::Constant(DIMENSION, STEP * (K - 1) / 2)) / ((STEP * (K - 1) / 2) * sqrt(DIMENSION)));
        EXPECT_TRUE(voctree::checkElements(features[j], "init"));
        featPtr.push_back(const_cast<FeatureFloat*> (&features.back()));
      }
    }

    initializer(featPtr, K, centers, voctree::L2<FeatureFloat,FeatureFloat>());

    // it's difficult to check the result as it is random, just check there are no weird things
    EXPECT_TRUE(voctree::checkVectorElements(centers, "initializer"));

  }
}

int main()
{
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
