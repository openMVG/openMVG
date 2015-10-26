#ifdef HAVE_CCTAG

#include "CCTagLocalizer.hpp"
#include "reconstructed_regions.hpp"
#include <openMVG/sfm/sfm_data_io.hpp>
#include <openMVG/matching/indMatch.hpp>
#include <openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp>

//@fixme move/redefine
#define POPART_COUT(x) std::cout << x << std::endl
#define POPART_CERR(x) std::cerr << x << std::endl

namespace openMVG {
namespace localization {

// inputs
// - sfmdata path
// - descriptorsFolder directory with the sift
// - vocTreeFilepath; 
// - weightsFilepath; 
bool CCTagLocalizer::init(const std::string &sfmFilePath,
                          const std::string &descriptorsFolder)
{
  using namespace openMVG::features;
  
  // load the sfm data containing the 3D reconstruction info
  if (!Load(_sfm_data, sfmFilePath, sfm::ESfM_Data::ALL)) 
  {
    std::cerr << std::endl
      << "The input SfM_Data file "<< sfmFilePath << " cannot be read." << std::endl;
    return false;
  }

  // this block is used to get the type of features (by default SIFT) used
  // for the reconstruction
  const std::string sImage_describer = stlplus::create_filespec(descriptorsFolder, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if(!regions_type)
  {
    POPART_CERR("Invalid: "
            << sImage_describer << " regions type file.");
    return false;
  }
  
  bool loadSuccessful = loadReconstructionDescriptors(_sfm_data, descriptorsFolder);
  
  if(!loadSuccessful)
  {
    POPART_CERR("Unable to load the descriptors");
    return false;
  }
  
  const sfm::Landmarks & landmarks = _sfm_data.GetLandmarks();
  
//  for(const auto & landmark : landmarks)
//  {
//    // Use the first observation to retrieve the associated descriptor.
//    const auto & firstObservation = *landmark.second.obs.begin();
//    
//    // Retrieve the Regions of the first observation
//    auto & reconstructedRegions = _regions_per_view[firstObservation.first];
//    
//    // Get the feature id: remap the index as we only load the reconstructed regions
//    const auto localFeatureId = reconstructedRegions._mapFullToLocal[firstObservation.second.id_feat];
//    
//    const auto & desc = reconstructedRegions._regions.Descriptors()[localFeatureId];
//    IndexT idCCTag = getCCTagId(desc);
//
//    // Insert <idCCTag, 3D point> into a map.
//    if (idCCTag!=UndefinedIndexT)
//    {
//      _cctagDatabase.emplace(idCCTag, landmark.second.X);
//    }
//  }
  return true;
}

/**
 * @brief Convert the descriptor representation into a CCTag ID.
 * @param[in] desc descriptor
 * @return cctag id or UndefinedIndexT if wrong cctag descriptor
 */
IndexT getCCTagId(const CCTagDescriptor & desc)
{
  std::size_t cctagId = UndefinedIndexT;
  for (int i = 0; i < desc.size(); ++i)
  {
    if (desc.getData()[i] == 1.0)
    {
      if (cctagId != UndefinedIndexT)
      {
        return UndefinedIndexT;
      }
      cctagId = i;
    }
    else if(desc.getData()[i] != 0)
    {
      return UndefinedIndexT;
    }
  }
  return cctagId;
}

bool CCTagLocalizer::loadReconstructionDescriptors(const sfm::SfM_Data & sfm_data,
                                                   const std::string & feat_directory)
{
  C_Progress_display my_progress_bar(sfm_data.GetViews().size(),
                                     std::cout, "\n- Regions Loading -\n");

  std::cout << "Build observations per view" << std::endl;
  // Build observations per view
  std::map<IndexT, std::vector<FeatureInImage> > observationsPerView;
  for(auto landmarkValue : sfm_data.structure)
  {
    IndexT trackId = landmarkValue.first;
    sfm::Landmark& landmark = landmarkValue.second;
    for(auto obs : landmark.obs)
    {
      const IndexT viewId = obs.first;
      const sfm::Observation& obs2d = obs.second;
      observationsPerView[viewId].push_back(FeatureInImage(obs2d.id_feat, trackId));
    }
  }
  for(auto featuresInImage : observationsPerView)
  {
    std::sort(featuresInImage.second.begin(), featuresInImage.second.end());
  }

  std::cout << "Load Features and Descriptors per view" << std::endl;
  // Read for each view the corresponding regions and store them
  for(sfm::Views::const_iterator iter = sfm_data.GetViews().begin();
          iter != sfm_data.GetViews().end(); ++iter, ++my_progress_bar)
  {
    const IndexT id_view = iter->second->id_view;
    Reconstructed_RegionsCCTag& reconstructedRegion = _regions_per_view[id_view];

    const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
    const std::string basename = stlplus::basename_part(sImageName);
    const std::string featFilepath = stlplus::create_filespec(feat_directory, basename, ".feat");
    const std::string descFilepath = stlplus::create_filespec(feat_directory, basename, ".desc");
    //    std::cout << "Feat: " << featFilepath << std::endl;
    //    std::cout << "Desc: " << descFilepath << std::endl;

    if(!reconstructedRegion._regions.Load(featFilepath, descFilepath))
    {
      std::cerr << "Invalid regions files for the view: " << sImageName << std::endl;
      return false;
    }

    // Filter descriptors to keep only the 3D reconstructed points
    reconstructedRegion.filterRegions(observationsPerView[id_view]);
  }
  return true;
}

bool CCTagLocalizer::localize(const image::Image<unsigned char> & imageGrey,
                const CCTagLocalizer::Parameters &param,
                bool useInputIntrinsics,
                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                LocalizationResult & localizationResult)
{
  // extract descriptors and features from image
  POPART_COUT("[features]\tExtract CCTag from query image");
  std::unique_ptr<features::Regions> tmpQueryRegions(new features::CCTAG_Regions());
  _image_describer.Describe(imageGrey, tmpQueryRegions);
  POPART_COUT("[features]\tExtract CCTAG done: found " << tmpQueryRegions->RegionCount() << " features");
  features::CCTAG_Regions queryRegions = *dynamic_cast<features::CCTAG_Regions*> (tmpQueryRegions.get());
  
  std::vector<IndexT> nearestKeyFrames;
  nearestKeyFrames.reserve(param._nNearestKeyFrames);
  
  kNearestKeyFrames(
          queryRegions,
          _regions_per_view,
          param._nNearestKeyFrames,
          nearestKeyFrames);
  
  // Set the minimum of the residual to infinite.
  double residualMin = std::numeric_limits<double>::max();
  IndexT indexBestKeyFrame = UndefinedIndexT;
  
  // Loop over all k nearest key frames in order to get the most geometrically 
  // consistent one.
  sfm::Image_Localizer_Match_Data bestResectionData;
  std::vector<pair<IndexT, IndexT> > bestAssociationIDs;
  geometry::Pose3 bestPose;
  
  for(const auto indexKeyFrame : nearestKeyFrames)
  {
    const Reconstructed_RegionsCCTag& matchedRegions = _regions_per_view[indexKeyFrame];
    
    // Matching
    std::vector<matching::IndMatch> vec_featureMatches;
    viewMatching(queryRegions, _regions_per_view[indexKeyFrame]._regions, vec_featureMatches);
    
    if ( vec_featureMatches.size() < 3 )
      continue;
    
    // D. recover the 2D-3D associations from the matches 
    // Each matched feature in the current similar image is associated to a 3D point,
    // hence we can recover the 2D-3D associations to estimate the pose
    // Prepare data for resection
    std::vector<pair<IndexT, IndexT> > associationIDsTemp;
    sfm::Image_Localizer_Match_Data resectionDataTemp;
    
    resectionDataTemp.error_max = param._errorMax;
    
    resectionDataTemp = sfm::Image_Localizer_Match_Data();
    resectionDataTemp.pt2D = Mat2X(2, vec_featureMatches.size());
    resectionDataTemp.pt3D = Mat3X(3, vec_featureMatches.size());
    
    // Get the 3D points associated to each matched feature
    std::size_t index = 0;
    for(const matching::IndMatch& featureMatch : vec_featureMatches)
    {
      assert(vec_featureMatches.size()>index);
      // the ID of the 3D point
      const IndexT trackId3D = matchedRegions._associated3dPoint[featureMatch._j];

      // prepare data for resectioning
      resectionDataTemp.pt3D.col(index) = _sfm_data.GetLandmarks().at(trackId3D).X;

      const Vec2 feat = queryRegions.GetRegionPosition(featureMatch._i);
      resectionDataTemp.pt2D.col(index) = feat;

      associationIDsTemp.emplace_back(trackId3D, featureMatch._i);
      ++index;
    }
    
    geometry::Pose3 poseTemp;
    
    bool bResection = sfm::SfM_Localizer::Localize(
            std::make_pair(imageGrey.Width(), imageGrey.Height()),
            // pass the input intrinsic if they are valid, null otherwise
            (useInputIntrinsics) ? &queryIntrinsics : nullptr,
            resectionDataTemp,
            poseTemp);

    if ( ( bResection ) && ( resectionDataTemp.error_max < residualMin) )
    {
      residualMin = resectionDataTemp.error_max;
      indexBestKeyFrame = indexKeyFrame;
      // Update best inital pose.
      bestPose = poseTemp;
      bestResectionData = resectionDataTemp;
      std::swap(bestAssociationIDs, associationIDsTemp);
    }
  }
  
  // If the resection has failed for all the nearest keyframes, the localization 
  // has failed.
  if ( indexBestKeyFrame == UndefinedIndexT ) 
  {
    return false;
  }
  
  // if we didn't use the provided intrinsics, estimate K from the projection
  // matrix estimated by the localizer and initialize the queryIntrinsics with
  // it and the image size. This will provide a first guess for the refine function
  if(!useInputIntrinsics)
  {
    // Decompose P matrix
    Mat3 K_, R_;
    Vec3 t_;
    // Decompose the projection matrix  to get K, R and t using 
    // RQ decomposition
    KRt_From_P(bestResectionData.projection_matrix, &K_, &R_, &t_);
    queryIntrinsics.setK(K_);
    queryIntrinsics.setWidth(imageGrey.Width());
    queryIntrinsics.setHeight(imageGrey.Height());
  }
  
  // Upper bound pixel(s) tolerance for residual errors
  double error_max = std::numeric_limits<double>::infinity();
  size_t max_iteration = 4096;
  
  // E. refine the estimated pose
  POPART_COUT("[poseEstimation]\tRefining estimated pose");
  bool refineStatus = sfm::SfM_Localizer::RefinePose(
          &queryIntrinsics, 
          bestPose, 
          bestResectionData, 
          true /*b_refine_pose*/, 
          param._refineIntrinsics /*b_refine_intrinsic*/);
  
  if(!refineStatus)
    POPART_COUT("[poseEstimation]\tRefine pose could not improve the estimation of the camera pose.");
  
  localizationResult = LocalizationResult(bestResectionData, bestAssociationIDs, bestPose, queryIntrinsics, true);

  return localizationResult.isValid();
  
 } 

CCTagLocalizer::CCTagLocalizer()
{
}

CCTagLocalizer::CCTagLocalizer(const CCTagLocalizer& orig)
{
}

CCTagLocalizer::~CCTagLocalizer()
{
}

void kNearestKeyFrames(
          const features::CCTAG_Regions & queryRegions,
          const CCTagRegionsPerViews & regionsPerView,
          std::size_t nNearestKeyFrames,
          std::vector<IndexT> & kNearestFrames)
{
  kNearestFrames.clear();
  
  // A std::multimap is used instead of a std::map because is very likely that the
  // similarity measure is equal for a subset of views in the CCTag regions case.
  std::multimap<float, IndexT> sortedViewSimilarities;
  
  for(const auto & keyFrame : regionsPerView)
  {
    const float similarity = viewSimilarity(queryRegions, keyFrame.second._regions);
    sortedViewSimilarities.emplace(similarity, keyFrame.first);
  }
  
  std::size_t counter = 0;
  for (auto rit = sortedViewSimilarities.crbegin(); rit != sortedViewSimilarities.crend(); ++rit)
  {
    kNearestFrames.push_back(rit->second);
    ++counter;
    
    if (counter == nNearestKeyFrames)
      break;
  }
}
 
void viewMatching(
        const features::CCTAG_Regions & regionsA,
        const features::CCTAG_Regions & regionsB,
        std::vector<matching::IndMatch> & vec_featureMatches)
{
  vec_featureMatches.clear();
  
  for(std::size_t i=0 ; i < regionsA.Descriptors().size() ; ++i)
  {
    const IndexT cctagIdA = getCCTagId(regionsA.Descriptors()[i]);
    // todo: Should be change to: Find in regionsB.Descriptors() the nearest 
    // descriptor to descriptorA. Currently, a cctag descriptor encode directly
    // the cctag id, then the id equality is tested.
    for(std::size_t j=0 ; j < regionsB.Descriptors().size() ; ++j)
    {
      const IndexT cctagIdB = getCCTagId(regionsB.Descriptors()[j]);
      if ( cctagIdA == cctagIdB )
      {
        vec_featureMatches.emplace_back(i,j);
        break;
      }
    }
  }
}
 
 
 
float viewSimilarity(
        const features::CCTAG_Regions & regionsA,
        const features::CCTAG_Regions & regionsB)
{
  assert(regionsA.DescriptorLength() == regionsB.DescriptorLength()); 
  
  const std::bitset<128> descriptorViewA = constructCCTagViewDescriptor(regionsA.Descriptors());
  const std::bitset<128> descriptorViewB = constructCCTagViewDescriptor(regionsB.Descriptors());
  
  // The similarity is the sum of all the cctags sharing the same id visible in both views.
  return (descriptorViewA & descriptorViewB).count();
}

std::bitset<128> constructCCTagViewDescriptor(
        const std::vector<CCTagDescriptor> & vCCTagDescriptors)
{
  std::bitset<128> descriptorView;
  for(const auto & cctagDescriptor : vCCTagDescriptors )
  {
    const IndexT cctagId = getCCTagId(cctagDescriptor);
    if ( cctagId != UndefinedIndexT)
    {
      descriptorView.set(cctagId, true);
    }
  }
  return descriptorView;
}

} // localization
} // openMVG

#endif //HAVE_CCTAG
