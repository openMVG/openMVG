
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "reconstructed_regions.hpp"

#include <openMVG/features/image_describer.hpp>
#include <nonFree/sift/SIFT_float_describer.hpp>
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/sfm/pipelines/localization/SfM_Localizer.hpp>
#include <openMVG/stl/stlMap.hpp>
#include <openMVG/voctree/vocabulary_tree.hpp>
#include <openMVG/voctree/database.hpp>
#include <openMVG/matching/matcher_kdtree_flann.hpp>
#include <openMVG/matching/regions_matcher.hpp>
#include <flann/algorithms/dist.h>


namespace openMVG {
namespace localization {

//@fixme find a better place or maje the class template?
typedef openMVG::features::Descriptor<float, 128> DescriptorFloat;
typedef Reconstructed_Regions<features::SIOPointFeature, float, 128> Reconstructed_RegionsT;

class VoctreeLocalizer
{

public:
  enum Algorithm : int { FirstBest=0, BestResult=1, AllResults=2, Cluster=3};
  static Algorithm initFromString(const std::string &value);
  
public:
  
  bool init(const std::string &sfmFilePath,
            const std::string &descriptorsFolder,
            const std::string &vocTreeFilepath,
            const std::string &weightsFilepath);
  
  // loadSfmData(const std::string & sfmDataPath)

  /**
   * @brief Load all the Descriptors who have contributed to the reconstruction.
   */
  bool loadReconstructionDescriptors(
    const sfm::SfM_Data & sfm_data,
    const std::string & feat_directory);
  
  /**
   * 
   * @param imageGray The input grayscale image
   * @param[in,out] queryIntrinsics Intrinsic parameters of the camera, they are used if the
   * flag useInputIntrinsics is set to true, otherwise they are estimated from the correspondences.
   * @param[in] numResults Number of images to retrieve from the vocabulary tree
   * @param[out] pose The camera pose
   * @param[in] useGuidedMatching Enable the guided match using the estimated fundamental
   * matrix 
   * @param[in] useInputIntrinsics Uses the \p queryIntrinsics as known calibration
   * @param[in] refineIntrinsics Refine the intrinsic parameters after resection
   * @param[in] algorithm Algorithm to use
   * @param[out] resection_data [optional] the 2D-3D correspondances used to compute the pose
   * @return true if the localization is successful
   */
  bool localize( const image::Image<unsigned char> & imageGray,
                cameras::Pinhole_Intrinsic &queryIntrinsics,
                const size_t numResults,
                geometry::Pose3 & pose,
                bool useGuidedMatching,
                bool useInputIntrinsics,
                bool refineIntrinsics,
                Algorithm algorithm,
                sfm::Image_Localizer_Match_Data * resection_data = nullptr);

  /**
   * @brief Try to localize an image in the database: it queries the database to 
   * retrieve \p numResults matching images and it tries to localize the query image
   * wrt the retrieve images in order of their score taking the first best result.
   *
   * @param imageGray The input grayscale image
   * @param[in,out] queryIntrinsics Intrinsic parameters of the camera, they are used if the
   * flag useInputIntrinsics is set to true, otherwise they are estimated from the correspondences.
   * @param[in] numResults Number of images to retrieve from the vocabulary tree
   * @param[out] pose The camera pose
   * @param[in] useGuidedMatching Enable the guided match using the estimated fundamental
   * matrix 
   * @param[in] useInputIntrinsics Uses the \p queryIntrinsics as known calibration
   * @param[in] refineIntrinsics Refine the intrinsic parameters after resection
   * @param[out] resection_data [optional] the 2D-3D correspondances used to compute the pose
   * @return true if the localization is successful
   */
  bool localizeFirstBestResult( const image::Image<unsigned char> & imageGray,
                cameras::Pinhole_Intrinsic &queryIntrinsics,
                const size_t numResults,
                geometry::Pose3 & pose,
                bool useGuidedMatching,
                bool useInputIntrinsics,
                bool refineIntrinsics,
                sfm::Image_Localizer_Match_Data * resection_data = nullptr);

  /**
   * @brief Try to localize an image in the database: it queries the database to 
   * retrieve \p numResults matching images and it tries to localize the query image
   * wrt the retrieve images in order of their score, collecting all the 2d-3d correspondences
   * and performing the resection with all these correspondences
   *
   * @param imageGray The input grayscale image
   * @param[in,out] queryIntrinsics Intrinsic parameters of the camera, they are used if the
   * flag useInputIntrinsics is set to true, otherwise they are estimated from the correspondences.
   * @param[in] numResults Number of images to retrieve from the vocabulary tree
   * @param[out] pose The camera pose
   * @param[in] useGuidedMatching Enable the guided match using the estimated fundamental
   * matrix 
   * @param[in] useInputIntrinsics Uses the \p queryIntrinsics as known calibration
   * @param[in] refineIntrinsics Refine the intrinsic parameters after resection
   * @param[out] resection_data [optional] the 2D-3D correspondences used to compute the pose
   * @return true if the localization is successful
   */
  bool localizeAllResults( const image::Image<unsigned char> & imageGray,
                cameras::Pinhole_Intrinsic &queryIntrinsics,
                const size_t numResults,
                geometry::Pose3 & pose,
                bool useGuidedMatching,
                bool useInputIntrinsics,
                bool refineIntrinsics,
                sfm::Image_Localizer_Match_Data * resection_data = nullptr);
  
  const sfm::SfM_Data& getSfMData() const {return _sfm_data; }

private:
  /**
   * @brief Load the vocabulary tree.
   */
  bool initDatabase(const std::string & vocTreeFilepath,
                                    const std::string & weightsFilepath,
                                    const std::string & feat_directory);

  typedef flann::L2<float> MetricT;
  typedef matching::ArrayMatcher_Kdtree_Flann<float, MetricT> MatcherT;
  bool robustMatching(matching::RegionsMatcherT<MatcherT> & matcher, 
                      const cameras::IntrinsicBase * queryIntrinsics,// the intrinsics of the image we are using as reference
                      const Reconstructed_RegionsT & regionsToMatch,
                      const cameras::IntrinsicBase * matchedIntrinsics,
                      const float fDistRatio,
                      const bool b_guided_matching,
                      const std::pair<size_t,size_t> & imageSizeI,     // size of the first image  
                      const std::pair<size_t,size_t> & imageSizeJ,     // size of the first image
                      std::vector<matching::IndMatch> & vec_featureMatches) const;
  
public:
  
  // for each view index, it contains the features and descriptors that have an
  // associated 3D point
  Hash_Map<IndexT, Reconstructed_RegionsT > _regions_per_view;
  
  // contains the 3D reconstruction data
  sfm::SfM_Data _sfm_data;
  
  // the feature extractor
  // @fixme do we want a generic image describer?
  features::SIFT_float_describer _image_describer;
  
  // the vocabulary tree used to generate the database and the visual images for
  // the query images
  voctree::VocabularyTree<DescriptorFloat> _voctree;
  
  // the database that stores the visual word representation of each image of
  // the original dataset
  voctree::Database _database;
  
  // this maps the docId in the database with the view index of the associated
  // image
  std::map<voctree::DocId, IndexT> _mapDocIdToView;
  
};

/**
 * @brief Print the name of the algorithm
 */
std::ostream& operator<<(std::ostream& os, VoctreeLocalizer::Algorithm a);

/**
 * @brief Get the type of algorithm from an integer
 */
std::istream& operator>>(std::istream &in, VoctreeLocalizer::Algorithm &a);


} // localization
} // openMVG
