#pragma once

#include "ILocalizer.hpp"
#include "LocalizationResult.hpp"
#include "VoctreeLocalizer.hpp"

#include <openMVG/features/features.hpp>
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
typedef Reconstructed_RegionsCCTag::FeatureT CCTagKeypoint;
typedef Hash_Map<IndexT, Reconstructed_RegionsCCTag > CCTagRegionsPerViews;

class CCTagLocalizer : public ILocalizer
{
  
  public:
  struct Parameters : LocalizerParameters
  {

    Parameters() : LocalizerParameters(), 
      _nNearestKeyFrames(4) { }
    
    size_t _nNearestKeyFrames;         //< number of best matching images to retrieve from the database                
  };
  
public:
  
  CCTagLocalizer(const std::string &sfmFilePath,
                 const std::string &descriptorsFolder);
   
 /**
   * @brief Just a wrapper around the different localization algorithm, the algorith
   * used to localized is chosen using \p param._algorithm
   * 
   * @param[in] imageGrey The input greyscale image.
   * @param[in] param The parameters for the localization.
   * @param[in] useInputIntrinsics Uses the \p queryIntrinsics as known calibration.
   * @param[in,out] queryIntrinsics Intrinsic parameters of the camera, they are used if the
   * flag useInputIntrinsics is set to true, otherwise they are estimated from the correspondences.
   * @param[out] localizationResult The localization result containing the pose and the associations.
   * @param[in] imagePath Optional complete path to the image, used only for debugging purposes.
   * @return  true if the image has been successfully localized.
   */
  bool localize(const image::Image<unsigned char> & imageGrey,
                const LocalizerParameters *parameters,
                bool useInputIntrinsics,
                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                LocalizationResult & localizationResult, const std::string& imagePath = std::string());

  bool localize(const std::unique_ptr<features::Regions> &queryRegions,
                const std::pair<std::size_t, std::size_t> &imageSize,
                const LocalizerParameters *parameters,
                bool useInputIntrinsics,
                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                LocalizationResult & localizationResult,
                const std::string& imagePath = std::string());
  /**
   * @brief Naive implementation of the localizer using the rig. Each image from
   * the rig is localized and then a bundle adjustment is run for optimizing the 
   * global pose.
   * 
   * @param[in] vec_imageGrey A vector containing all the images from the rig
   * @param[in] parameters The parameters for the localization.
   * @param[in,out] vec_queryIntrinsics Vector containing the intrinsic parameters of the cameras
   * @param[in] vec_subPoses A vector containing the N-1 subposes of each camera wrt the main camera
   * @param[out] rigPose The pose of the rig expressed as the pose of the main camera
   * @return true if the rig has been successfully localized.
   */
  bool localizeRig(const std::vector<image::Image<unsigned char> > & vec_imageGrey,
                   const LocalizerParameters *parameters,
                   std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                   const std::vector<geometry::Pose3 > &vec_subPoses,
                   geometry::Pose3 &rigPose);
  

  bool localizeRig(const std::vector<std::unique_ptr<features::Regions> > & vec_queryRegions,
                   const std::vector<std::pair<std::size_t, std::size_t> > &imageSize,
                   const LocalizerParameters *param,
                   std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                   const std::vector<geometry::Pose3 > &vec_subPoses,
                   geometry::Pose3 &rigPose);
  
  
  /**
   * @brief Given the input Regions, it retrieves all the 2D-3D associations from
   * the nearest k-frames in the database. The associations are retrived in terms
   * of region index and 3D point index along with the number of times (\p occurrences) that the 
   * pair has been found. \p pt2D and \p pt3D contains the coordinates of the corresponding
   * points of the associations, in the same order as in \p occurences.
   * 
   * @param[in] queryRegions The input query regions containing the extracted 
   * markers from the query image
   * @param[in] param The parameters to use
   * @param[out] occurences A map with a pair of indices for each association as 
   * key and its occurrence as value
   * @param[out] pt2D The set of 2D points of the associations as they are given in \p queryRegions
   * @param[out] pt3D The set of 3D points of the associations
   */
  void getAllAssociations(const features::CCTAG_Regions &queryRegions,
                          const CCTagLocalizer::Parameters &param,
                          std::map< std::pair<IndexT, IndexT>, std::size_t > &occurences,
                          Mat &pt2D,
                          Mat &pt3D) const;
  
  virtual ~CCTagLocalizer();

private:
  
  bool loadReconstructionDescriptors(
    const sfm::SfM_Data & sfm_data,
    const std::string & feat_directory);
  
  // for each view index, it contains the cctag features and descriptors that have an
  // associated 3D point
  CCTagRegionsPerViews _regions_per_view;
   
  // the feature extractor
  features::CCTAG_Image_describer _image_describer;
  
  //
  //std::map<IndexT, Vec3> _cctagDatabase;
};

 /**
   * @brief Retrieve the k nearest views in a collection of views based on a query
   *        consisting in a set of CCTag regions.
   * 
   * @param[in] queryRegions Set of CCTag regions in the query.
   * @param[in] Collection of views containing a set of cctag regions.
   * @param[in] nNearestKeyFrames Number of nearest neighbours to return.
   * @param[out] kNearestFrames Set of computed indices associated to the k nearest views.
   * @param[in] similarityThreshold A threshold to retrieve only the kframes having 
  *  at least \p similarityThreshold similarity score.
   */
void kNearestKeyFrames(
          const features::CCTAG_Regions & queryRegions,
          const CCTagRegionsPerViews & _regions_per_view,
          std::size_t nNearestKeyFrames,
          std::vector<IndexT> & kNearestFrames,
          const float similarityThreshold = .0f);
/**
 * @brief Given a set of CCTag descriptors seen in a view, it creates a descriptor for the view: the
 * view descriptor is a 128 bit array (ie the number of possible markers) whose 
 * bits are 0 or 1 whether the corresponding marker is seen or not. E.g. if the 
 * bit in position 8 is 1 it means that the marker with ID 8 has been seen by the view.
 * 
 * @param[in] vCCTagDescriptors The input descriptors associated to the view.
 * @return The view descriptor as a set of bit representing the visibility of
 * each possible marker for that view.
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

