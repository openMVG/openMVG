#ifdef HAVE_CCTAG

#pragma once

#include <openMVG/localization/VoctreeLocalizer.hpp>

#include "openMVG/features/features.hpp"
#include "LocalizationResult.hpp"
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
                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                LocalizationResult & localizationResult);
  
  bool localize(const std::vector<image::Image<unsigned char> > & vec_imageGrey,
                const Parameters &param,
                const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                const std::vector<geometry::Pose3 > &vec_subPoses,
                geometry::Pose3 rigPose);
  
  bool localize(const std::vector<std::unique_ptr<features::Regions> > & vec_queryRegions,
                const Parameters &param,
                const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                const std::vector<geometry::Pose3 > &vec_subPoses,
                geometry::Pose3 rigPose);
  
  CCTagLocalizer();

  virtual ~CCTagLocalizer();
  
  const sfm::SfM_Data& getSfMData() const {return _sfm_data; }
  
  void getAllAssociationsFromNearestKFrames(const features::CCTAG_Regions &queryRegions,
                                            const CCTagLocalizer::Parameters &param,
                                            std::map< pair<IndexT, IndexT>, pair<Vec3, Vec2> > &associations) const;
  
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
/**
 * @brief Given a set of CCTag descriptors seen in a view, it creates a descriptor for the view: the
 * view descriptor is a 128 bit array (ie the number of possible markers) whose 
 * bits are 0 or 1 whether the corresponding marker is seen or not. E.g. if the 
 * bit in position 8 is 1 it means that the marker with ID 8 has been seen by the view
 * 
 * @param[in] vCCTagDescriptors The input descriptors associated to the view
 * @return The view descriptor as a set of bit representing the visibility of 
 * each possible marker for that view
 */
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

#endif //HAVE_CCTAG