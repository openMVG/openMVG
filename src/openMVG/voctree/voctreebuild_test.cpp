/* 
 * File:   voctreebuild_test.cpp
 * Author: sgaspari
 *
 * Created on February 11, 2016, 2:33 PM
 */

#include <openMVG/voctree/tree_builder.hpp>

#include <testing/testing.h>

#include <Eigen/Core>

#include <iostream>
#include <fstream>
#include <vector>


TEST(voctree, voctreeBuilder)
{
  using namespace openMVG;

  const std::string treeName = "test.tree";

  const std::size_t DIMENSION = 3;
  const std::size_t FEATURENUMBER = 100;

  const float kepsf = 10e-8;
  const std::size_t K = 4;
  const std::size_t LEVELS = 3;
  const std::size_t LEAVESNUMBER = std::pow(K, LEVELS);


  const std::size_t STEP = 1;

  typedef Eigen::Matrix<float, 1, DIMENSION> FeatureFloat;
  typedef std::vector<FeatureFloat, Eigen::aligned_allocator<FeatureFloat> > FeatureFloatVector;
  typedef std::vector<FeatureFloat* > FeaturePointerVector;

  // generate a random vector of features
  FeatureFloatVector features;
  features.reserve(FEATURENUMBER * LEAVESNUMBER);

  for(std::size_t i = 0; i < LEAVESNUMBER; ++i)
  {
    // at each i iteration translate the cluster by STEP*i
    for(std::size_t j = 0; j < FEATURENUMBER; ++j)
    {
      features.push_back((FeatureFloat::Random(1, DIMENSION) + Eigen::MatrixXf::Constant(1, DIMENSION, STEP * i) - Eigen::MatrixXf::Constant(1, DIMENSION, STEP * (LEAVESNUMBER - 1) / 2)) / ((STEP * (LEAVESNUMBER - 1) / 2) * sqrt(DIMENSION)));
      EXPECT_TRUE(voctree::checkElements(features[j], "init"));
    }
  }

  // build the tree
  voctree::TreeBuilder<FeatureFloat> builder(FeatureFloat::Zero());
  builder.setVerbose(0);
  builder.kmeans().setRestarts(10);
  std::cout << "Building a tree of L = " << LEVELS << " levels with a branching factor of k = " << K << std::endl;
  builder.build(features, K, LEVELS);
  std::cout << builder.tree().centers().size() << " centers" << std::endl;

//  // if we quantize the features, considering each created cluster as a document, 
//  // every feature in the cluster should be assigned to the same visual word
//  for(std::size_t i = 0; i < LEAVESNUMBER; ++i)
//  {
//    // get the visual word of the first element
//    const voctree::Word &vword = builder.tree().quantize(features[i * FEATURENUMBER ]);
//    for(std::size_t j = 1; j < FEATURENUMBER; ++j)
//    {
//      // it should be true but it is not always the case
////      BOOST_WARN_NE(vword, builder.tree().quantize(features[i * FEATURENUMBER + j ]));
//      if(vword != builder.tree().quantize(features[i * FEATURENUMBER + j ]))
//      {
//        std::cerr << "Warning: feature " << j << " not assigned to his visual word " 
//                << i << " (visual word "<< vword << " instead)" << std::endl; 
//      }
//    }
//  }

  // the centers should all be valid in this configuration
  std::vector<uint8_t> valid = builder.tree().validCenters();
  for(std::size_t i = 0; i < valid.size(); ++i)
    EXPECT_TRUE(valid[i] != 0);

  builder.tree().save(treeName);

  voctree::MutableVocabularyTree<FeatureFloat> loadedtree;
  loadedtree.load(treeName);

  // check the centers are the same
  FeatureFloatVector centerOrig = builder.tree().centers();
  FeatureFloatVector centerLoad = loadedtree.centers();

  EXPECT_EQ(centerOrig.size(), centerLoad.size());

  voctree::L2<FeatureFloat, FeatureFloat> distance;
  for(std::size_t i = 0; i < centerOrig.size(); ++i)
  {
    EXPECT_NEAR(distance(centerOrig[i], centerLoad[i]), 0, kepsf);
  }


//  voctree::printFeatVector( features ); 
}

int main()
{
  TestResult tr;
  return TestRegistry::runAllTests(tr);
}

