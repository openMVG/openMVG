
// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "reconstructed_regions.hpp"
#include "LocalizationResult.hpp"
#include "ILocalizer.hpp"

#include <openMVG/features/image_describer.hpp>
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

class VoctreeLocalizer : public ILocalizer
{

public:
  enum Algorithm : int { FirstBest=0, BestResult=1, AllResults=2, Cluster=3};
  static Algorithm initFromString(const std::string &value);
  
public:
  struct Parameters : LocalizerParameters
  {

    Parameters() : LocalizerParameters(), 
      _useGuidedMatching(false),
      _algorithm(Algorithm::FirstBest),
      _numResults(4),
      _maxResults(10),
      _numCommonViews(3) { }
    
    bool _useGuidedMatching;    //< Enable/disable guided matching when matching images
    Algorithm _algorithm;       //< algorithm to use for localization
    size_t _numResults;         //< number of best matching images to retrieve from the database
    size_t _maxResults;         
    size_t _numCommonViews;     //< number minimum common images in which a point must be seen to be used in cluster tracking
  };
  
public:
  
  /**
   * @brief Initialize a localizer based on a vocabulary tree
   * 
   * @param[in] sfmFilePath The path to the sfmdata file containing the scene 
   * reconstruction.
   * @param[in] descriptorsFolder The path to the directory containing the features 
   * of the scene (.desc and .feat files).
   * @param[in] vocTreeFilepath The path to the vocabulary tree (usually a .tree file).
   * @param[in] weightsFilepath Optional path to the weights of the vocabulary 
   * tree (usually a .weights file), if not provided the weights will be recomputed 
   * when all the documents are added.
   * @param[in] useSIFT_CCTAG Optional and enabled only if the CCTAG are available. 
   * It enable the use of combined SIFT and CCTAG features.
   */
  VoctreeLocalizer(const std::string &sfmFilePath,
                   const std::string &descriptorsFolder,
                   const std::string &vocTreeFilepath,
                   const std::string &weightsFilepath
#ifdef HAVE_CCTAG
                   , bool useSIFT_CCTAG
#endif
                  );
  
  /**
   * @brief Just a wrapper around the different localization algorithm, the algorith
   * used to localized is chosen using \p param._algorithm
   * 
   * @param[in] imageGray The input greyscale image
   * @param[in] param The parameters for the localization
   * @param[in] useInputIntrinsics Uses the \p queryIntrinsics as known calibration
   * @param[in,out] queryIntrinsics Intrinsic parameters of the camera, they are used if the
   * flag useInputIntrinsics is set to true, otherwise they are estimated from the correspondences.
   * @param[out] pose The camera pose
   * @param[out] resection_data the 2D-3D correspondences used to compute the pose
   * @return true if the localization is successful
   */
  bool localize(const image::Image<unsigned char> & imageGrey,
                const LocalizerParameters *param,
                bool useInputIntrinsics,
                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                LocalizationResult &localizationResult);

  /**
   * @brief Try to localize an image in the database: it queries the database to 
   * retrieve \p numResults matching images and it tries to localize the query image
   * wrt the retrieve images in order of their score taking the first best result.
   *
   * @param[in] imageGray The input greayscale image
   * @param[in] param The parameters for the localization
   * @param[in] useInputIntrinsics Uses the \p queryIntrinsics as known calibration
   * @param[in,out] queryIntrinsics Intrinsic parameters of the camera, they are used if the
   * flag useInputIntrinsics is set to true, otherwise they are estimated from the correspondences.
   * @param[out] pose The camera pose
   * @param[out] resection_data the 2D-3D correspondences used to compute the pose
   * @param[out] associationIDs the ids of the 2D-3D correspondences used to compute the pose
   * @return true if the localization is successful
   */
  bool localizeFirstBestResult(const image::Image<unsigned char> & imageGrey,
                const Parameters &param,
                bool useInputIntrinsics,
                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                LocalizationResult &localizationResult);

  /**
   * @brief Try to localize an image in the database: it queries the database to 
   * retrieve \p numResults matching images and it tries to localize the query image
   * wrt the retrieve images in order of their score, collecting all the 2d-3d correspondences
   * and performing the resection with all these correspondences
   *
   * @param[in] imageGray The input greyscale image
   * @param[in] param The parameters for the localization
   * @param[in] useInputIntrinsics Uses the \p queryIntrinsics as known calibration
   * @param[in,out] queryIntrinsics Intrinsic parameters of the camera, they are used if the
   * flag useInputIntrinsics is set to true, otherwise they are estimated from the correspondences.
   * @param[out] pose The camera pose
   * @param[out] resection_data the 2D-3D correspondences used to compute the pose
   * @param[out] associationIDs the ids of the 2D-3D correspondences used to compute the pose
   * @return true if the localization is successful
   */
  bool localizeAllResults(const image::Image<unsigned char> & imageGrey,
                const Parameters &param,
                bool useInputIntrinsics,
                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                LocalizationResult &localizationResult);

private:
  /**
   * @brief Load the vocabulary tree.

   * @param[in] vocTreeFilepath The path to the directory containing the features 
   * of the scene (.desc and .feat files).
   * @param[in] weightsFilepath weightsFilepath Optional path to the weights of the vocabulary 
   * tree (usually a .weights file), if not provided the weights will be recomputed 
   * when all the documents are added.
   * @param[in] feat_directory The path to the directory containing the features 
   * of the scene (.desc and .feat files).
   * @return true if everything went ok
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
  
  /**
   * @brief Load all the Descriptors who have contributed to the reconstruction.
   * deprecated.. now inside initDatabase
   */
  bool loadReconstructionDescriptors(
    const sfm::SfM_Data & sfm_data,
    const std::string & feat_directory);
  
  
public:
  
  // for each view index, it contains the features and descriptors that have an
  // associated 3D point
  Hash_Map<IndexT, Reconstructed_RegionsT > _regions_per_view;
  
  // the feature extractor
  // @fixme do we want a generic image describer?
//  features::SIFT_float_describer _image_describer;
  features::Image_describer* _image_describer;
  
  // the vocabulary tree used to generate the database and the visual images for
  // the query images
  voctree::VocabularyTree<DescriptorFloat> _voctree;
  
  // the database that stores the visual word representation of each image of
  // the original dataset
  voctree::Database _database;
  
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
