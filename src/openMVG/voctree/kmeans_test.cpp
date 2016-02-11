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

int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
