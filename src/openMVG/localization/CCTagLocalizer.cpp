#ifdef HAVE_CCTAG

#include "CCTagLocalizer.hpp"
#include "reconstructed_regions.hpp"
#include <openMVG/sfm/sfm_data_io.hpp>
#include <openMVG/matching/indMatch.hpp>
#include <openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp>

//@fixme move/redefine
#define POPART_COUT(x) std::cout << x << std::endl
#define POPART_CERR(x) std::cerr << x << std::endl
#define POPART_COUT_DEBUG(x) std::cout << x << std::endl

namespace openMVG {
namespace localization {

CCTagLocalizer::CCTagLocalizer(const std::string &sfmFilePath,
                               const std::string &descriptorsFolder)
{
  using namespace openMVG::features;
  
  // load the sfm data containing the 3D reconstruction info
  if (!Load(_sfm_data, sfmFilePath, sfm::ESfM_Data::ALL)) 
  {
    std::cerr << std::endl
      << "The input SfM_Data file "<< sfmFilePath << " cannot be read." << std::endl;
    throw std::invalid_argument("The input SfM_Data file "+ sfmFilePath + " cannot be read.");
  }

  // this block is used to get the type of features (by default SIFT) used
  // for the reconstruction
  const std::string sImage_describer = stlplus::create_filespec(descriptorsFolder, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if(!regions_type)
  {
    POPART_CERR("Invalid: "
            << sImage_describer << " regions type file.");
    throw std::invalid_argument("Invalid: "+ sImage_describer + " regions type file.");
  }
  
  bool loadSuccessful = loadReconstructionDescriptors(_sfm_data, descriptorsFolder);
  
  if(!loadSuccessful)
  {
    POPART_CERR("Unable to load the descriptors");
    throw std::invalid_argument("Unable to load the descriptors from "+descriptorsFolder);
  }
  
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
  
  _isInit = true;
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
  std::vector<bool> presentIds(128,false); // @todo Assume a maximum library size of 128 unique ids.
  // Read for each view the corresponding regions and store them
  for(sfm::Views::const_iterator iter = sfm_data.GetViews().begin();
          iter != sfm_data.GetViews().end(); ++iter, ++my_progress_bar)
  {
    const IndexT id_view = iter->second->id_view;
    Reconstructed_RegionsCCTag& reconstructedRegion = _regions_per_view[id_view];

    const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
    std::string featFilepath = stlplus::create_filespec(feat_directory, std::to_string(iter->first), ".feat");
    std::string descFilepath = stlplus::create_filespec(feat_directory, std::to_string(iter->first), ".desc");

    if(!(stlplus::is_file(featFilepath) && stlplus::is_file(descFilepath)))
    {
      // legacy compatibility, if the features are not named using the UID convention
      // let's try with the old-fashion naming convention
      const std::string basename = stlplus::basename_part(sImageName);
      featFilepath = stlplus::create_filespec(feat_directory, basename, ".feat");
      descFilepath = stlplus::create_filespec(feat_directory, basename, ".desc");
      if(!(stlplus::is_file(featFilepath) && stlplus::is_file(descFilepath)))
      {
        POPART_CERR("Cannot find the features for image " << sImageName 
                << " neither using the UID naming convention nor the image name based convention");
        return false;
      }
    }

    if(!reconstructedRegion._regions.Load(featFilepath, descFilepath))
    {
      std::cerr << "Invalid regions files for the view: " << sImageName << std::endl;
      return false;
    }

    // Filter descriptors to keep only the 3D reconstructed points
    reconstructedRegion.filterCCTagRegions(observationsPerView[id_view]);
    
    // Update the visibility mask
    reconstructedRegion.updateLandmarksVisibility(presentIds);
  }

  // Display the cctag ids over all cctag landmarks present in the database
  std::cout << std::endl << "Present CCTag landmarks present in the database: " << std::endl;
  for(std::size_t i = 0; i < presentIds.size(); ++i)
  {
    if (presentIds[i])
      std::cout <<  i+1 << " ";
  }

  std::cout <<  std::endl << std::endl;
  
  return true;
}

bool CCTagLocalizer::localize(const image::Image<unsigned char> & imageGrey,
                const LocalizerParameters *parameters,
                bool useInputIntrinsics,
                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                LocalizationResult & localizationResult)
{
  const CCTagLocalizer::Parameters *param = static_cast<const CCTagLocalizer::Parameters *>(parameters);
  if(!param)
  {
    // error!
    throw std::invalid_argument("The parameters are not in the right format!!");
  }
  // extract descriptors and features from image
  POPART_COUT("[features]\tExtract CCTag from query image");
  std::unique_ptr<features::Regions> tmpQueryRegions(new features::CCTAG_Regions());
  _image_describer.Describe(imageGrey, tmpQueryRegions);
  POPART_COUT("[features]\tExtract CCTAG done: found " << tmpQueryRegions->RegionCount() << " features");
  features::CCTAG_Regions &queryRegions = *dynamic_cast<features::CCTAG_Regions*> (tmpQueryRegions.get());
  
  std::vector<IndexT> nearestKeyFrames;
  nearestKeyFrames.reserve(param->_nNearestKeyFrames);
  
  kNearestKeyFrames(
          queryRegions,
          _regions_per_view,
          param->_nNearestKeyFrames,
          nearestKeyFrames);
  
  // Set the minimum of the residual to infinite.
  double residualMin = std::numeric_limits<double>::max();
  IndexT indexBestKeyFrame = UndefinedIndexT;
  
  // Loop over all k nearest key frames in order to get the most geometrically 
  // consistent one.
  sfm::Image_Localizer_Match_Data bestResectionData;
  std::vector<pair<IndexT, IndexT> > bestAssociationIDs;
  geometry::Pose3 bestPose;
  
  POPART_COUT_DEBUG("nearestKeyFrames.size() = " << nearestKeyFrames.size());
  for(const IndexT indexKeyFrame : nearestKeyFrames)
  {
    POPART_COUT_DEBUG("[localization]\tProcessing nearest kframe " << indexKeyFrame);
    POPART_COUT_DEBUG(_sfm_data.GetViews().at(indexKeyFrame)->s_Img_path);
    const Reconstructed_RegionsCCTag& matchedRegions = _regions_per_view[indexKeyFrame];
    
    // Matching
    std::vector<matching::IndMatch> vec_featureMatches;
    viewMatching(queryRegions, _regions_per_view[indexKeyFrame]._regions, vec_featureMatches);
    
    if ( vec_featureMatches.size() < 3 )
    {
      POPART_COUT("[localization]\tSkipping kframe " << indexKeyFrame << " as it contains only "<< vec_featureMatches.size()<<" matches");
      continue;
    }
    
    // D. recover the 2D-3D associations from the matches 
    // Each matched feature in the current similar image is associated to a 3D point,
    // hence we can recover the 2D-3D associations to estimate the pose
    // Prepare data for resection
    std::vector<pair<IndexT, IndexT> > associationIDsTemp;
    sfm::Image_Localizer_Match_Data resectionDataTemp;
    
    resectionDataTemp.error_max = param->_errorMax;
    
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
  bestResectionData.error_max = std::numeric_limits<double>::infinity();
  bestResectionData.max_iteration = 4096;
  
  // E. refine the estimated pose
  POPART_COUT("[poseEstimation]\tRefining estimated pose");
  bool refineStatus = sfm::SfM_Localizer::RefinePose(
          &queryIntrinsics, 
          bestPose, 
          bestResectionData, 
          true /*b_refine_pose*/, 
          param->_refineIntrinsics /*b_refine_intrinsic*/);
  
  if(!refineStatus)
    POPART_COUT("[poseEstimation]\tRefine pose could not improve the estimation of the camera pose.");
  
  localizationResult = LocalizationResult(bestResectionData, bestAssociationIDs, bestPose, queryIntrinsics, true);

  return localizationResult.isValid();
  
 } 

CCTagLocalizer::~CCTagLocalizer()
{
}


bool CCTagLocalizer::localize(const std::vector<image::Image<unsigned char> > & vec_imageGrey,
                              const Parameters &param,
                              const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                              const std::vector<geometry::Pose3 > &vec_subPoses,
                              geometry::Pose3 rigPose)
{
  assert(vec_imageGrey.size()==vec_queryIntrinsics.size());
  assert(vec_imageGrey.size()==vec_subPoses.size());

  const size_t numCams = vec_imageGrey.size();

  std::vector<std::unique_ptr<features::Regions> > vec_queryRegions;
  vec_queryRegions.reserve(numCams);

  for( size_t i = 0; i < numCams; ++i )
  {
    // initialize the regions
    vec_queryRegions.emplace_back(new features::CCTAG_Regions());
  }

  //@todo parallelize?
  for( size_t i = 0; i < numCams; ++i )
  {
    // for each image, extract descriptors and features 
    POPART_COUT("[features]\tExtract CCTag from query image " << i);
    _image_describer.Describe(vec_imageGrey[i], vec_queryRegions[i]);
    POPART_COUT("[features]\tExtract CCTAG done: found " << vec_queryRegions[i]->RegionCount() << " features for query image " << i);
  }

  return localize(vec_queryRegions, param, vec_queryIntrinsics, vec_subPoses, rigPose);
}


bool CCTagLocalizer::localize(const std::vector<std::unique_ptr<features::Regions> > & vec_queryRegions,
              const Parameters &param,
              const std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
              const std::vector<geometry::Pose3 > &vec_subPoses,
              geometry::Pose3 rigPose)
{
  assert(vec_queryRegions.size()==vec_queryIntrinsics.size());
  assert(vec_queryRegions.size()==vec_subPoses.size());   

  const size_t numCams = vec_queryRegions.size();

  vector<std::map< pair<IndexT, IndexT>, pair<Vec3, Vec2> > > vec_associations(numCams);

  // for each camera retrieve the associations
  //@todo parallelize?
  size_t numAssociations = 0;
  for( size_t i = 0; i < numCams; ++i )
  {

    // this map is used to collect the 2d-3d associations as we go through the images
    // the key is a pair <Id3D, Id2d>
    // the element is the pair 3D point - 2D point
    std::map< pair<IndexT, IndexT>, pair<Vec3, Vec2> > &associations = vec_associations[i];
    features::CCTAG_Regions &queryRegions = *dynamic_cast<features::CCTAG_Regions*> (vec_queryRegions[i].get());
    getAllAssociationsFromNearestKFrames(queryRegions, param, associations);
    numAssociations += associations.size();
  }

  const size_t minNumAssociations = 5;  //possible parameter?
  if(numAssociations < minNumAssociations)
  {
    POPART_COUT("[poseEstimation]\tonly " << numAssociations << " have been found, not enough to do the resection!");
    return false;
  }

  // prepare resection data
  sfm::Image_Localizer_Match_Data resectionData;
  resectionData.pt2D = Mat2X(2, numAssociations);
  resectionData.pt3D = Mat3X(3, numAssociations);

  // fill the point matrices with the relevant points
  size_t index = 0;
  for( size_t i = 0; i < numCams; ++i )
  {
    std::map< pair<IndexT, IndexT>, pair<Vec3, Vec2> > &associations = vec_associations[i];

    for(const auto &ass : associations)
    {
      const cameras::Pinhole_Intrinsic_Radial_K3 &currCamera = vec_queryIntrinsics[i];
      const geometry::Pose3 &currPose = vec_subPoses[i];
       // recopy all the points in the matching structure
      //@todo 
      // undistort and normalize the 2D points 
      // we first remove the distortion and then we transform the undistorted point in
      // normalized camera coordinates (inv(K)*undistortedPoint)
      resectionData.pt2D.col(index) = currCamera.ima2cam(currCamera.remove_disto(ass.second.second));
      // multiply the 3D point by the camera rototranslation
      resectionData.pt3D.col(index) = currPose(ass.second.first);

      ++index;
    }

  }
  assert(index==numAssociations);

  // do the resection with all the associations
  // the resection in this case must be done in normalized coordinates for the 
  // 2D feature points (ie inv(K)*feat) so that all the 2D points are independent
  // from their camera parameters (f, pp, and distortion)
  // estimate the pose
  // Do the resectioning: compute the camera pose.
  resectionData.error_max = param._errorMax;
  
  // this is a 'fake' intrinsics for the camera rig: all the 2D points are already
  // undistorted and in their normalized camera coordinates, hence the K is the 
  // identity matrix (f=1, pp=(0,0)), image size is irrelevant/not used so set to 0
  cameras::Pinhole_Intrinsic rigIntrinsics(0,0,1,0,0);
  POPART_COUT("[poseEstimation]\tEstimating camera pose...");
  bool bResection = sfm::SfM_Localizer::Localize(std::make_pair(0,0), // image size is not used for calibrated case
                                                 &rigIntrinsics,
                                                 resectionData,
                                                 rigPose);

  if(!bResection)
  {
    POPART_COUT("[poseEstimation]\tResection FAILED");
    return false;
  }
  POPART_COUT("[poseEstimation]\tResection SUCCEDED");

  POPART_COUT("R est\n" << rigPose.rotation());
  POPART_COUT("t est\n" << rigPose.translation());
  
  // E. refine the estimated pose
  POPART_COUT("[poseEstimation]\tRefining estimated pose");
  bool refineStatus = sfm::SfM_Localizer::RefinePose(&rigIntrinsics, 
                                                     rigPose, 
                                                     resectionData, 
                                                     true /*b_refine_pose*/, 
                                                     false /*b_refine_intrinsic*/);
  if(!refineStatus)
    POPART_COUT("Refine pose could not improve the estimation of the camera pose.");

  {
    // just temporary code to evaluate the estimated pose @todo remove it
    POPART_COUT("R refined\n" << rigPose.rotation());
    POPART_COUT("t refined\n" << rigPose.translation());
    POPART_COUT("K refined\n" << rigIntrinsics.K());
  }
  
  //@fixme return something meaningful
  return true;
}


void CCTagLocalizer::getAllAssociationsFromNearestKFrames(const features::CCTAG_Regions &queryRegions,
                                                          const CCTagLocalizer::Parameters &param,
                                                          std::map< pair<IndexT, IndexT>, pair<Vec3, Vec2> > &associations) const
{
  std::vector<IndexT> nearestKeyFrames;
  nearestKeyFrames.reserve(param._nNearestKeyFrames);
  
  kNearestKeyFrames(
          queryRegions,
          _regions_per_view,
          param._nNearestKeyFrames,
          nearestKeyFrames);
  
  POPART_COUT_DEBUG("nearestKeyFrames.size() = " << nearestKeyFrames.size());
  for(const IndexT indexKeyFrame : nearestKeyFrames)
  {
    POPART_COUT_DEBUG(indexKeyFrame);
    POPART_COUT_DEBUG(_sfm_data.GetViews().at(indexKeyFrame)->s_Img_path);
    const Reconstructed_RegionsCCTag& matchedRegions = _regions_per_view.at(indexKeyFrame);
    
    // Matching
    std::vector<matching::IndMatch> vec_featureMatches;
    viewMatching(queryRegions, _regions_per_view.at(indexKeyFrame)._regions, vec_featureMatches);
    
    // D. recover the 2D-3D associations from the matches 
    // Each matched feature in the current similar image is associated to a 3D point,
    // hence we can recover the 2D-3D associations to estimate the pose
    
    // Get the 3D points associated to each matched feature
    std::size_t index = 0;
    for(const matching::IndMatch& featureMatch : vec_featureMatches)
    {
      assert(vec_featureMatches.size()>index);
      // the ID of the 3D point
      const IndexT pt3D_id = matchedRegions._associated3dPoint[featureMatch._j];
      const IndexT pt2D_id = featureMatch._i;
      
      const auto key = std::make_pair(pt3D_id, pt2D_id);
      if(associations.count(key))
      {
        // we already have this association, skip it
        continue;
      }
      
      // prepare data for resectioning
      const auto &point3d = _sfm_data.GetLandmarks().at(pt3D_id).X;

      const Vec2 feat = queryRegions.GetRegionPosition(pt2D_id);
      
      associations.insert(std::make_pair(key, std::make_pair(point3d, feat)));
      
      ++index;
    }
  }
  
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
    const IndexT cctagIdA = features::getCCTagId(regionsA.Descriptors()[i]);
    // todo: Should be change to: Find in regionsB.Descriptors() the nearest 
    // descriptor to descriptorA. Currently, a cctag descriptor encode directly
    // the cctag id, then the id equality is tested.
    for(std::size_t j=0 ; j < regionsB.Descriptors().size() ; ++j)
    {
      const IndexT cctagIdB = features::getCCTagId(regionsB.Descriptors()[j]);
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
    const IndexT cctagId = features::getCCTagId(cctagDescriptor);
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
