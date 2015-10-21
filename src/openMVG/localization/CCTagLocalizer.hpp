#pragma once

#include <openMVG/localization/localization.hpp>

#include "openMVG/features/features.hpp"
#include <openMVG/features/image_describer.hpp>
#include <openMVG/features/cctag/CCTAG_describer.hpp>
#include <openMVG/sfm/sfm_data.hpp>
#include <openMVG/sfm/pipelines/localization/SfM_Localizer.hpp>

#include <iostream>
#include <bitset>

namespace openMVG {
namespace localization {

typedef Reconstructed_Regions<features::SIOPointFeature, unsigned char, 128> Reconstructed_RegionsCCTag; 
typedef Reconstructed_RegionsCCTag::DescriptorT CCTagDescriptor;
typedef Hash_Map<IndexT, Reconstructed_RegionsCCTag > CCTagRegionsPerViews;

class CCTagLocalizer {
  
  public:
  struct Parameters 
  {

    Parameters() :
      _useGuidedMatching(false),
      _refineIntrinsics(false),
      _nNearestKeyFrames(4),
      _numCommonViews(3),
      _fDistRatio(0.6),
      _featurePreset(features::EDESCRIBER_PRESET::ULTRA_PRESET),
      _errorMax(std::numeric_limits<double>::max()) { }
    
    bool _useGuidedMatching;    //< Enable/disable guided matching when matching images
    bool _refineIntrinsics;     //< whether or not the Intrinsics of the query camera has to be refined
    size_t _nNearestKeyFrames;         //< number of best matching images to retrieve from the database        
    size_t _numCommonViews;     //< number minimum common images in which a point must be seen to be used in cluster tracking
    float _fDistRatio;          //< the ratio distance to use when matching feature with the ratio test
    features::EDESCRIBER_PRESET _featurePreset; //< the preset to use for feature extraction of the query image
    double _errorMax;           
  };
  
public:
  
  bool init(const std::string &sfmFilePath,
            const std::string &descriptorsFolder);
  
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
                const Parameters &param,
                bool useInputIntrinsics,
                cameras::Pinhole_Intrinsic &queryIntrinsics,
                geometry::Pose3 & pose,
                sfm::Image_Localizer_Match_Data &resection_data,
                std::vector<pair<IndexT, IndexT> > &associationIDs);
  
  
  CCTagLocalizer();
  CCTagLocalizer(const CCTagLocalizer& orig);
  virtual ~CCTagLocalizer();
  const sfm::SfM_Data& getSfMData() const {return _sfm_data; }
private:
  
  bool loadReconstructionDescriptors(
    const sfm::SfM_Data & sfm_data,
    const std::string & feat_directory);
  
  // for each view index, it contains the cctag features and descriptors that have an
  // associated 3D point
  
  CCTagRegionsPerViews _regions_per_view;
  // contains the 3D reconstruction data
  sfm::SfM_Data _sfm_data;
  
  // the feature extractor
  features::CCTAG_Image_describer _image_describer;
  
  //
  //std::map<IndexT, Vec3> _cctagDatabase;
};

IndexT getCCTagId(const CCTagDescriptor & desc);

 /**
   * @brief Retrieve the k nearest views in a collection of views based on a query
   *        consisting in a set of CCTag regions.
   * 
   * @param[in] queryRegions Set of CCTag regions in the query
   * @param[in] Collection of views containing a set of cctag regions 
   * @param[in] nNearestKeyFrames Number of nearest neighbours to return
   * @param[out] kNearestFrames Set of computed indices associated to the k nearest views
   */
void kNearestKeyFrames(
          const features::CCTAG_Regions & queryRegions,
          const CCTagRegionsPerViews & _regions_per_view,
          std::size_t nNearestKeyFrames,
          std::vector<IndexT> & kNearestFrames);

std::bitset<128> constructCCTagViewDescriptor(
        const std::vector<CCTagDescriptor> & vCCTagDescriptors);

float viewSimilarity(
        const features::CCTAG_Regions & regionsA,
        const features::CCTAG_Regions & regionsB);

void viewMatching(
        const features::CCTAG_Regions & regionsA,
        const features::CCTAG_Regions & regionsB,
        std::vector<matching::IndMatch> & vec_featureMatches);

} // namespace localization
} // openMVG
