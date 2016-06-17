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

  const std::size_t DIMENSION = 128;
  const std::size_t FEATURENUMBER = 500;

  const std::size_t K = 10;

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

  // now try to generate k cluster well far away and compare
  features.clear();
  featPtr.clear();
  features.reserve(FEATURENUMBER * K);
  featPtr.reserve(features.size());
  for(std::size_t i = 0; i < K; ++i)
  {
    // at each i iteration translate the cluster by 5*i
    for(std::size_t j = 0; j < FEATURENUMBER; ++j)
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

  const int FEATURENUMBER = 500;
  const std::size_t numTrial = 3;
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
TEST(kmeans, kmeanSimple)
{
  using namespace openMVG;

  std::cout << "Testing kmeans..." << std::endl;

  const std::size_t DIMENSION = 8;
  const std::size_t FEATURENUMBER = 500;

  const std::size_t K = 30;

  const std::size_t STEP = 5 * K;

  typedef Eigen::Matrix<float, 1, DIMENSION> FeatureFloat;
  typedef std::vector<FeatureFloat, Eigen::aligned_allocator<FeatureFloat> > FeatureFloatVector;
  typedef std::vector<FeatureFloat* > FeaturePointerVector;

  // generate a random vector of features
  FeatureFloatVector features;
  FeatureFloatVector centers;
  FeatureFloatVector centersGT;
  std::vector<unsigned int> membership;
  features.reserve(FEATURENUMBER);
  centers.reserve(K);

  voctree::SimpleKmeans<FeatureFloat> kmeans(FeatureFloat::Zero());
  kmeans.setVerbose(0);
  kmeans.setRestarts(5);

  for(std::size_t trial = 0; trial < 10; ++trial)
  {
    // now try to generate k cluster well far away and comapare
    features.clear();
    membership.clear();
    centersGT.clear();
    centers.clear();
    features.reserve(FEATURENUMBER * K);
    membership.reserve(features.size());
    centersGT.reserve(K);
    centers.reserve(K);
    for(std::size_t i = 0; i < K; ++i)
    {
      // at each i iteration translate the cluster by STEP*i
      for(std::size_t j = 0; j < FEATURENUMBER; ++j)
      {
        features.push_back((FeatureFloat::Random(1, DIMENSION) + Eigen::MatrixXf::Constant(1, DIMENSION, STEP * i) - Eigen::MatrixXf::Constant(1, DIMENSION, STEP * (K - 1) / 2)) / ((STEP * (K - 1) / 2) * sqrt(DIMENSION)));
        EXPECT_TRUE(voctree::checkElements(features[j], "init"));
      }
      centersGT.push_back((Eigen::MatrixXf::Constant(1, DIMENSION, STEP * i) - Eigen::MatrixXf::Constant(1, DIMENSION, STEP * (K - 1) / 2)) / ((STEP * (K - 1) / 2) * sqrt(DIMENSION)));
    }

    voctree::SimpleKmeans<FeatureFloat>::squared_distance_type dist = kmeans.cluster(features, K, centers, membership);

//    voctree::printFeatVector( centers );

    voctree::L2<FeatureFloat, FeatureFloat> distance;
    for(std::size_t i = 0; i < K; ++i)
    {
      voctree::SimpleKmeans<FeatureFloat>::squared_distance_type bestDist = std::numeric_limits<voctree::SimpleKmeans<FeatureFloat>::squared_distance_type>::max();
      for(std::size_t j = 0; j < K; ++j)
      {
        voctree::SimpleKmeans<FeatureFloat>::squared_distance_type centerDist = distance(centers[j], centersGT[i]);
        if(centerDist < bestDist)
          bestDist = centerDist;
      }
    }

    std::vector<unsigned int> h(K, 0);
    for(std::size_t i = 0; i < membership.size(); ++i)
    {
      ++h[membership[i]];
    }
    for(std::size_t i = 0; i < h.size(); ++i)
    {
      EXPECT_TRUE(h[i] > 0);
    }
  }
}

TEST(kmeans, kmeanVarying)
{
  using namespace openMVG;
  std::cout << "Testing kmeans with variable k and DIM..." << std::endl;

  const std::size_t FEATURENUMBER = 300;
  const std::size_t numTrial = 3;

  // generate random values for K and DIMENSION
  std::default_random_engine generator;
  std::uniform_int_distribution<std::size_t> dimGen(3, 128);
  std::uniform_int_distribution<std::size_t> kGen(6, 300);

  for(std::size_t trial = 0; trial < numTrial; ++trial)
  {
    const std::size_t DIMENSION = dimGen(generator);
    const std::size_t K = kGen(generator);
    const std::size_t STEP = 5 * K;
    std::cout << "\tTrial " << trial + 1 << "/" << numTrial << " with K = " << K << " and DIMENSION = " << DIMENSION << std::endl;

    typedef Eigen::RowVectorXf FeatureFloat;
    typedef std::vector<FeatureFloat, Eigen::aligned_allocator<FeatureFloat> > FeatureFloatVector;
    typedef std::vector<FeatureFloat* > FeaturePointerVector;

    // generate a random vector of features
    FeatureFloatVector features;
    FeatureFloatVector centers;
    FeatureFloatVector centersGT;
    std::vector<unsigned int> membership;

    voctree::SimpleKmeans<FeatureFloat> kmeans(FeatureFloat::Zero(DIMENSION));
    kmeans.setVerbose(0);
    kmeans.setRestarts(3);

    features.reserve(FEATURENUMBER * K);
    membership.reserve(features.size());
    centersGT.reserve(K);
    centers.reserve(K);

    // now try to generate k cluster well far away and comapare
    for(std::size_t i = 0; i < K; ++i)
    {
      // at each i iteration translate the cluster by STEP*i
      for(std::size_t j = 0; j < FEATURENUMBER; ++j)
      {
        features.push_back((FeatureFloat::Random(DIMENSION) + FeatureFloat::Constant(DIMENSION, STEP * i) - FeatureFloat::Constant(DIMENSION, STEP * (K - 1) / 2)) / ((STEP * (K - 1) / 2) * sqrt(DIMENSION)));
        EXPECT_TRUE(voctree::checkElements(features[j], "init"));
      }
      centersGT.push_back((FeatureFloat::Constant(DIMENSION, STEP * i) - FeatureFloat::Constant(DIMENSION, STEP * (K - 1) / 2)) / ((STEP * (K - 1) / 2) * sqrt(DIMENSION)));
    }

    voctree::SimpleKmeans<FeatureFloat>::squared_distance_type dist = kmeans.cluster(features, K, centers, membership);

//    voctree::printFeatVector( features );
//    voctree::printFeatVector(centers);
//    voctree::printFeatVector(centersGT);

    voctree::L2<FeatureFloat, FeatureFloat> distance;
    voctree::SimpleKmeans<FeatureFloat>::squared_distance_type globDist = 0.0;
    for(size_t i = 0; i < K; ++i)
    {
      voctree::SimpleKmeans<FeatureFloat>::squared_distance_type bestDist = std::numeric_limits<voctree::SimpleKmeans<FeatureFloat>::squared_distance_type>::max();
      for(std::size_t j = 0; j < K; ++j)
      {
        voctree::SimpleKmeans<FeatureFloat>::squared_distance_type centerDist = distance(centers[j], centersGT[i]);
        if(centerDist < bestDist)
          bestDist = centerDist;
      }
      globDist += bestDist;
    }
    std::cout << "center distance " << globDist << std::endl;

    std::vector<size_t> h(K, 0);
    for(size_t i = 0; i < membership.size(); ++i)
    {
      ++h[membership[i]];
    }
    for(size_t i = 0; i < h.size(); ++i)
    {
//      std::cout << h[i] << std::endl;
      EXPECT_TRUE(h[i] > 0);
      EXPECT_EQ(h[i], FEATURENUMBER);
    }
  }

}

int main()
{
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}
