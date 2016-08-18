#include "CCTagLocalizer.hpp"
#include "reconstructed_regions.hpp"
#include "optimization.hpp"
#include "rigResection.hpp"

#include <openMVG/features/svgVisualization.hpp>
#include <openMVG/sfm/sfm_data_io.hpp>
#include <openMVG/matching/indMatch.hpp>
#include <openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp>
#include <openMVG/system/timer.hpp>
#include <openMVG/logger.hpp>

#include <cctag/ICCTag.hpp>

#include <boost/filesystem.hpp>

#include <algorithm>

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
    POPART_CERR("\n\nIf the error says \"JSON Parsing failed - provided NVP not found\" "
        "it's likely that you have to convert your sfm_data to a recent version supporting "
        "polymorphic Views. You can run the python script convertSfmData.py to update an existing sfmdata.");
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
  std::vector<int> counterCCtagsInImage = {0, 0, 0, 0, 0, 0};
  // Read for each view the corresponding regions and store them
  for(const auto &iter : sfm_data.GetViews())
  {
    const IndexT id_view = iter.second->id_view;
    Reconstructed_RegionsCCTag& reconstructedRegion = _regions_per_view[id_view];

    const std::string &sImageName = iter.second.get()->s_Img_path;
    std::string featFilepath = stlplus::create_filespec(feat_directory, std::to_string(iter.first), ".feat");
    std::string descFilepath = stlplus::create_filespec(feat_directory, std::to_string(iter.first), ".desc");

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
    
    ++my_progress_bar;
  }

  {
    // just debugging stuff -- print for each image the visible reconstructed cctag
    // and create an histogram of cctags per image
    for(const auto &iter : sfm_data.GetViews())
    {
      const IndexT id_view = iter.second->id_view;
      Reconstructed_RegionsCCTag& reconstructedRegion = _regions_per_view[id_view];
      const std::string &sImageName = iter.second.get()->s_Img_path;

      std::cout << "Image " << sImageName;
      if(reconstructedRegion._regions.Descriptors().size() == 0 )
      {
        counterCCtagsInImage[0] +=1;
        std::cout << " does not contain any cctag!!!";
      }
      else
      {
        std::cout << " contains CCTag Id:\t";
        for(const auto &desc : reconstructedRegion._regions.Descriptors())
        {
          const IndexT cctagIdA = features::getCCTagId(desc);
          if(cctagIdA != UndefinedIndexT)
            std::cout << cctagIdA << " ";
        }
        // Update histogram
        int countcctag = reconstructedRegion._regions.Descriptors().size();
        if(countcctag >= 5)
          counterCCtagsInImage[5] +=1;
        else
          counterCCtagsInImage[countcctag] += 1;
      }
      std::cout << "\n";
    }
    
    // Display histogram
    std::cout << std::endl << "Histogram of number of cctags in images :" << std::endl;
    for(std::size_t i = 0; i < 5; i++)
      std::cout << "Images with " << i << "  CCTags : " << counterCCtagsInImage[i] << std::endl;
    std::cout << "Images with 5+ CCTags : " << counterCCtagsInImage[5] << std::endl << std::endl;

    // Display the cctag ids over all cctag landmarks present in the database
    std::cout << std::endl << "Found " << std::count(presentIds.begin(), presentIds.end(), true) 
            << " CCTag in the database with " << _sfm_data.GetLandmarks().size() << " associated 3D points\n"
            "The CCTag id in the database are: " << std::endl;
    for(std::size_t i = 0; i < presentIds.size(); ++i)
    {
      if(presentIds[i])
        std::cout << i + 1 << " ";
    }

    std::cout << std::endl << std::endl;

  }
  
  return true;
}

bool CCTagLocalizer::localize(const image::Image<unsigned char> & imageGrey,
                              const LocalizerParameters *parameters,
                              bool useInputIntrinsics,
                              cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                              LocalizationResult & localizationResult, 
                              const std::string& imagePath)
{
  namespace bfs = boost::filesystem;
  
  const CCTagLocalizer::Parameters *param = static_cast<const CCTagLocalizer::Parameters *>(parameters);
  if(!param)
  {
    throw std::invalid_argument("The CCTag localizer parameters are not in the right format.");
  }
  // extract descriptors and features from image
  POPART_COUT("[features]\tExtract CCTag from query image");
  std::unique_ptr<features::Regions> tmpQueryRegions(new features::CCTAG_Regions());
  _image_describer.Set_configuration_preset(param->_featurePreset);
  _image_describer.Describe(imageGrey, tmpQueryRegions);
  POPART_COUT("[features]\tExtract CCTAG done: found " << tmpQueryRegions->RegionCount() << " features");
  
  std::pair<std::size_t, std::size_t> imageSize = std::make_pair(imageGrey.Width(),imageGrey.Height());
  
  if(!param->_visualDebug.empty() && !imagePath.empty())
  {
    // it automatically throws an exception if the cast does not work
    features::CCTAG_Regions &queryRegions = *dynamic_cast<features::CCTAG_Regions*> (tmpQueryRegions.get());
    
    // just debugging -- save the svg image with detected cctag
    features::saveCCTag2SVG(imagePath, 
                            imageSize, 
                            queryRegions, 
                            param->_visualDebug+"/"+bfs::path(imagePath).stem().string()+".svg");
  }

//  return localize(tmpQueryRegions,
  return localizeAllAssociations(tmpQueryRegions,
                  imageSize,
                  parameters,
                  useInputIntrinsics,
                  queryIntrinsics,
                  localizationResult,
                  imagePath);
}

bool CCTagLocalizer::localizeAllAssociations(const std::unique_ptr<features::Regions> &genQueryRegions,
                                              const std::pair<std::size_t, std::size_t> &imageSize,
                                              const LocalizerParameters *parameters,
                                              bool useInputIntrinsics,
                                              cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                                              LocalizationResult & localizationResult,
                                              const std::string& imagePath)
{
  namespace bfs = boost::filesystem;
  
  const CCTagLocalizer::Parameters *param = static_cast<const CCTagLocalizer::Parameters *>(parameters);
  if(!param)
  {
    throw std::invalid_argument("The CCTag localizer parameters are not in the right format.");
  }
  
  // it automatically throws an exception if the cast does not work
  features::CCTAG_Regions &queryRegions = *dynamic_cast<features::CCTAG_Regions*> (genQueryRegions.get());
  
  std::map< std::pair<IndexT, IndexT>, std::size_t > occurences;
  sfm::Image_Localizer_Match_Data resectionData;
  
  system::Timer timer;
  getAllAssociations(queryRegions, imageSize, *param, occurences, resectionData.pt2D, resectionData.pt3D, imagePath);
  POPART_COUT("[Matching]\tRetrieving associations took " << timer.elapsedMs() << "ms");
  
  const std::size_t numCollectedPts = occurences.size();
  
  // create an vector of <feat3D_id, feat2D_id>
  std::vector<pair<IndexT, IndexT> > associationIDs;
  associationIDs.reserve(numCollectedPts);

  for(const auto &ass : occurences)
  {
    // recopy the associations IDs in the vector
    associationIDs.push_back(ass.first);
  }
  
  assert(associationIDs.size() == numCollectedPts);
  assert(resectionData.pt2D.cols() == numCollectedPts);
  assert(resectionData.pt3D.cols() == numCollectedPts);

  geometry::Pose3 pose;
  
  timer.reset();
  // estimate the pose
  resectionData.error_max = param->_errorMax;
  POPART_COUT("[poseEstimation]\tEstimating camera pose...");
  const bool bResection = sfm::SfM_Localizer::Localize(imageSize,
                                                      // pass the input intrinsic if they are valid, null otherwise
                                                      (useInputIntrinsics) ? &queryIntrinsics : nullptr,
                                                      resectionData,
                                                      pose,
                                                      param->_resectionEstimator);
  
  if(!bResection)
  {
    POPART_COUT("[poseEstimation]\tResection failed");
    if(!param->_visualDebug.empty() && !imagePath.empty())
    {
//      namespace bfs = boost::filesystem;
//      features::saveFeatures2SVG(imagePath,
//                                 imageSize,
//                                 resectionData.pt2D,
//                                 param._visualDebug + "/" + bfs::path(imagePath).stem().string() + ".associations.svg");
    }
    localizationResult = LocalizationResult();
    return localizationResult.isValid();
  }
  POPART_COUT("[poseEstimation]\tResection SUCCEDED");

  POPART_COUT("R est\n" << pose.rotation());
  POPART_COUT("t est\n" << pose.translation());

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
    KRt_From_P(resectionData.projection_matrix, &K_, &R_, &t_);
    queryIntrinsics.setK(K_);
    POPART_COUT("K estimated\n" << K_);
    queryIntrinsics.setWidth(imageSize.first);
    queryIntrinsics.setHeight(imageSize.second);
  }

  // refine the estimated pose
  POPART_COUT("[poseEstimation]\tRefining estimated pose");
  const bool b_refine_pose = true;
  const bool refineStatus = sfm::SfM_Localizer::RefinePose(&queryIntrinsics,
                                                            pose,
                                                            resectionData,
                                                            b_refine_pose,
                                                            param->_refineIntrinsics);
  if(!refineStatus)
    POPART_COUT("Refine pose failed.");

  if(!param->_visualDebug.empty() && !imagePath.empty())
  {
//    namespace bfs = boost::filesystem;
//    features::saveFeatures2SVG(imagePath,
//                              imageSize,
//                              resectionData.pt2D,
//                              param._visualDebug + "/" + bfs::path(imagePath).stem().string() + ".associations.svg",
//                              &resectionData.vec_inliers);
  }
  
  POPART_COUT("[poseEstimation]\tPose estimation took " << timer.elapsedMs() << "ms.");

  //@todo for now just empty, to be added to getAllAssociations
  std::vector<voctree::DocMatch> matchedImages;
  localizationResult = LocalizationResult(resectionData, associationIDs, pose, queryIntrinsics, matchedImages, refineStatus);

  {
    // just debugging this block can be safely removed or commented out
    POPART_COUT("R refined\n" << pose.rotation());
    POPART_COUT("t refined\n" << pose.translation());
    POPART_COUT("K refined\n" << queryIntrinsics.K());

    const Mat2X residuals = localizationResult.computeInliersResiduals();

    const auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();
    POPART_COUT("RMSE = " << std::sqrt(sqrErrors.mean())
                << " min = " << std::sqrt(sqrErrors.minCoeff())
                << " max = " << std::sqrt(sqrErrors.maxCoeff()));
  }

  return localizationResult.isValid();
  
  
}


bool CCTagLocalizer::localize(const std::unique_ptr<features::Regions> &genQueryRegions,
                              const std::pair<std::size_t, std::size_t> &imageSize,
                              const LocalizerParameters *parameters,
                              bool useInputIntrinsics,
                              cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                              LocalizationResult & localizationResult,
                              const std::string& imagePath)
{
  namespace bfs = boost::filesystem;
  
  const CCTagLocalizer::Parameters *param = static_cast<const CCTagLocalizer::Parameters *>(parameters);
  if(!param)
  {
    throw std::invalid_argument("The CCTag localizer parameters are not in the right format.");
  }
  
  // it automatically throws an exception if the cast does not work
  features::CCTAG_Regions &queryRegions = *dynamic_cast<features::CCTAG_Regions*> (genQueryRegions.get());
   
  std::vector<IndexT> nearestKeyFrames;
  nearestKeyFrames.reserve(param->_nNearestKeyFrames);
  
  kNearestKeyFrames(queryRegions,
                    _regions_per_view,
                    param->_nNearestKeyFrames,
                    nearestKeyFrames, 4);
  
  // Set the minimum of the residual to infinite.
  double residualMin = std::numeric_limits<double>::max();
  IndexT indexBestKeyFrame = UndefinedIndexT;
  
  // Loop over all k nearest key frames in order to get the most geometrically 
  // consistent one.
  sfm::Image_Localizer_Match_Data bestResectionData;
    // Upper bound pixel(s) tolerance for residual errors
  bestResectionData.error_max = std::numeric_limits<double>::infinity();
  bestResectionData.max_iteration = 4096;
  
  std::vector<pair<IndexT, IndexT> > bestAssociationIDs;
  geometry::Pose3 bestPose;
  
  POPART_COUT_DEBUG("nearestKeyFrames.size() = " << nearestKeyFrames.size());
  for(const IndexT indexKeyFrame : nearestKeyFrames)
  {
    POPART_COUT("[matching]\tProcessing nearest kframe " << indexKeyFrame 
            << " (" << _sfm_data.GetViews().at(indexKeyFrame)->s_Img_path << ")");
    const Reconstructed_RegionsCCTag& matchedRegions = _regions_per_view[indexKeyFrame];
    
    // Matching
    std::vector<matching::IndMatch> vec_featureMatches;
    viewMatching(queryRegions, _regions_per_view[indexKeyFrame]._regions, vec_featureMatches);
    
    if(!param->_visualDebug.empty() && !imagePath.empty())
    {
      const sfm::View *mview = _sfm_data.GetViews().at(indexKeyFrame).get();
      const std::string queryImage = bfs::path(imagePath).stem().string();
      const std::string matchedImage = bfs::path(mview->s_Img_path).stem().string();
      const std::string matchedPath = (bfs::path(_sfm_data.s_root_path) /  bfs::path(mview->s_Img_path)).string();

      // the directory where to save the feature matches
      const auto baseDir = bfs::path(param->_visualDebug) / queryImage;
      if((!bfs::exists(baseDir)))
      {
        POPART_COUT("created " << baseDir.string());
        bfs::create_directories(baseDir);
      }
      
      // the final filename for the output svg file as a composition of the query
      // image and the matched image
      auto outputName = baseDir / queryImage;
      outputName += "_";
      outputName += matchedImage;
      outputName += ".svg";
      
      const bool showNotMatched = true;
      features::saveCCTagMatches2SVG(imagePath, 
                                     imageSize, 
                                     queryRegions,
                                     matchedPath,
                                     std::make_pair(mview->ui_width, mview->ui_height), 
                                     _regions_per_view[indexKeyFrame]._regions,
                                     vec_featureMatches,
                                     outputName.string(),
                                     showNotMatched ); 
    }
    
    if ( vec_featureMatches.size() < 3 )
    {
      POPART_COUT("[matching]\tSkipping kframe " << indexKeyFrame << " as it contains only "<< vec_featureMatches.size()<<" matches");
      continue;
    }
    POPART_COUT("[matching]\tFound "<< vec_featureMatches.size()<<" matches");
    
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
    
    bool bResection = sfm::SfM_Localizer::Localize(imageSize,
                                                  // pass the input intrinsic if they are valid, null otherwise
                                                  (useInputIntrinsics) ? &queryIntrinsics : nullptr,
                                                  resectionDataTemp,
                                                  poseTemp,
                                                  param->_resectionEstimator);

    if ( ( bResection ) && ( resectionDataTemp.error_max < residualMin) )
    {
      residualMin = resectionDataTemp.error_max;
      indexBestKeyFrame = indexKeyFrame;
      // Update best initial pose.
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
    queryIntrinsics.setWidth(imageSize.first);
    queryIntrinsics.setHeight(imageSize.second);
  }
  
  // E. refine the estimated pose
  POPART_COUT("[poseEstimation]\tRefining estimated pose");
  const bool b_refine_pose = true;
  bool refineStatus = sfm::SfM_Localizer::RefinePose(&queryIntrinsics, 
                                                     bestPose, 
                                                     bestResectionData, 
                                                     b_refine_pose, 
                                                     param->_refineIntrinsics);
  
  if(!refineStatus)
    POPART_COUT("[poseEstimation]\tRefine pose failed.");

  localizationResult = LocalizationResult(bestResectionData, bestAssociationIDs, bestPose, queryIntrinsics, std::vector<voctree::DocMatch>(), refineStatus);

  return localizationResult.isValid();
  
 } 

CCTagLocalizer::~CCTagLocalizer()
{
}

// subposes is n-1 as we consider the first camera as the main camera and the 
// reference frame of the grid
bool CCTagLocalizer::localizeRig(const std::vector<image::Image<unsigned char> > & vec_imageGrey,
                                 const LocalizerParameters *parameters,
                                 std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                                 const std::vector<geometry::Pose3 > &vec_subPoses,
                                 geometry::Pose3 &rigPose,
                                 std::vector<localization::LocalizationResult> & vec_locResults)
{
  const CCTagLocalizer::Parameters *param = static_cast<const CCTagLocalizer::Parameters *>(parameters);
  if(!param)
  {
    throw std::invalid_argument("The CCTag localizer parameters are not in the right format.");
  }
  const size_t numCams = vec_imageGrey.size();
  assert(numCams == vec_queryIntrinsics.size());
  assert(numCams == vec_subPoses.size() - 1);

  std::vector<std::unique_ptr<features::Regions> > vec_queryRegions(numCams);
  std::vector<std::pair<std::size_t, std::size_t> > vec_imageSize;
  
  //@todo parallelize?
  for(size_t i = 0; i < numCams; ++i)
  {
    // extract descriptors and features from each image
    vec_queryRegions[i] = std::unique_ptr<features::Regions>(new features::CCTAG_Regions());
    POPART_COUT("[features]\tExtract CCTag from query image...");
    _image_describer.Set_configuration_preset(param->_featurePreset);
    _image_describer.Describe(vec_imageGrey[i], vec_queryRegions[i]);
    POPART_COUT("[features]\tExtract CCTAG done: found " <<  vec_queryRegions[i]->RegionCount() << " features");
    // add the image size for this image
    vec_imageSize.emplace_back(vec_imageGrey[i].Width(), vec_imageGrey[i].Height());
  }
  assert(vec_imageSize.size() == vec_queryRegions.size());
          
  return localizeRig(vec_queryRegions,
                     vec_imageSize,
                     parameters,
                     vec_queryIntrinsics,
                     vec_subPoses,
                     rigPose,
                     vec_locResults);
}

bool CCTagLocalizer::localizeRig(const std::vector<std::unique_ptr<features::Regions> > & vec_queryRegions,
                                 const std::vector<std::pair<std::size_t, std::size_t> > &vec_imageSize,
                                 const LocalizerParameters *parameters,
                                 std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                                 const std::vector<geometry::Pose3 > &vec_subPoses,
                                 geometry::Pose3 &rigPose,
                                 std::vector<LocalizationResult>& vec_locResults)
{
#ifdef HAVE_OPENGV
  if(!parameters->_useLocalizeRigNaive)
  {
    POPART_COUT("Using localizeRig_naive()");
    return localizeRig_opengv(vec_queryRegions,
                              vec_imageSize,
                              parameters,
                              vec_queryIntrinsics,
                              vec_subPoses,
                              rigPose,
                              vec_locResults);
  }
  else
#endif
  {
    if(!parameters->_useLocalizeRigNaive)
      POPART_COUT("OpenGV is not available. Fallback to localizeRig_naive().");
    return localizeRig_naive(vec_queryRegions,
                             vec_imageSize,
                             parameters,
                             vec_queryIntrinsics,
                             vec_subPoses,
                             rigPose,
                             vec_locResults);
  }
}

#ifdef HAVE_OPENGV
bool CCTagLocalizer::localizeRig_opengv(const std::vector<std::unique_ptr<features::Regions> > & vec_queryRegions,
                                 const std::vector<std::pair<std::size_t, std::size_t> > &imageSize,
                                 const LocalizerParameters *parameters,
                                 std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                                 const std::vector<geometry::Pose3 > &vec_subPoses,
                                 geometry::Pose3 &rigPose,
                                 std::vector<LocalizationResult>& vec_locResults)
{
  const size_t numCams = vec_queryRegions.size();
  assert(numCams == vec_queryIntrinsics.size());
  assert(numCams == vec_subPoses.size() - 1);   

  const CCTagLocalizer::Parameters *param = static_cast<const CCTagLocalizer::Parameters *>(parameters);
  if(!param)
  {
    throw std::invalid_argument("The CCTag localizer parameters are not in the right format.");
  }

  vector<std::map< pair<IndexT, IndexT>, std::size_t > > vec_occurrences(numCams);
  vector<Mat> vec_pts3D(numCams);
  vector<Mat> vec_pts2D(numCams);

  // for each camera retrieve the associations
  //@todo parallelize?
  size_t numAssociations = 0;
  for( size_t i = 0; i < numCams; ++i )
  {

    // this map is used to collect the 2d-3d associations as we go through the images
    // the key is a pair <Id3D, Id2d>
    // the element is the pair 3D point - 2D point
    auto &occurrences = vec_occurrences[i];
    Mat &pts3D = vec_pts3D[i];
    Mat &pts2D = vec_pts2D[i];
    features::CCTAG_Regions &queryRegions = *dynamic_cast<features::CCTAG_Regions*> (vec_queryRegions[i].get());
    getAllAssociations(queryRegions, imageSize[i],*param, occurrences, pts2D, pts3D);
    numAssociations += occurrences.size();
  }
  
  // @todo Here it could be possible to filter the associations according to their
  // occurrences, eg giving priority to those associations that are more frequent

  const size_t minNumAssociations = 5;  //possible parameter?
  if(numAssociations < minNumAssociations)
  {
    POPART_COUT("[poseEstimation]\tonly " << numAssociations << " have been found, not enough to do the resection!");
    return false;
  }

  std::vector<std::vector<std::size_t> > vec_inliers;
  const bool resectionOk = rigResection(vec_pts2D,
                                        vec_pts3D,
                                        vec_queryIntrinsics,
                                        vec_subPoses,
                                        rigPose, 
                                        vec_inliers);

  if(!resectionOk)
  {
    for(std::size_t cam = 0; cam < numCams; ++cam)
    {
      // empty result with isValid set to false
      vec_locResults.emplace_back();
    }
    return resectionOk;
  }
  
  const bool refineOk = refineRigPose(vec_pts2D,
                                      vec_pts3D,
                                      vec_inliers,
                                      vec_queryIntrinsics,
                                      vec_subPoses,
                                      rigPose);
  
  vec_locResults.clear();
  vec_locResults.reserve(numCams);
  
  // create localization results
  for(std::size_t cam = 0; cam < numCams; ++cam)
  {

    const auto &intrinsics = vec_queryIntrinsics[cam];

    // compute the (absolute) pose of each camera: for the main camera it's the 
    // rig pose, for the others, combine the subpose with the rig pose
    geometry::Pose3 pose;
    if(cam == 0)
    {
      pose = rigPose;
    }
    else
    {
      // main camera: q1 ~ [R1 t1] Q = [I 0] A   where A = [R1 t1] Q  
      // another camera: q2 ~ [R2 t2] Q = [R2 t2]*inv([R1 t1]) A 
      // and subPose12 = [R12 t12] = [R2 t2]*inv([R1 t1])
      // With rigResection() we compute [R1 t1] (aka rigPose), hence:
      // subPose12 = [R12 t12] = [R2 t2]*inv([R1 t1]) and we need [R2 t2], ie the absolute pose
      // => [R1 t1] * subPose12 = [R2 t2]
      // => rigPose * subPose12 = [R2 t2]
      pose = vec_subPoses[cam] * rigPose;
    }
    
    // create matchData
    sfm::Image_Localizer_Match_Data matchData;
    matchData.vec_inliers = vec_inliers[cam];
    matchData.error_max = param->_errorMax;
    matchData.projection_matrix = intrinsics.get_projective_equivalent(pose);
    matchData.pt2D = vec_pts2D[cam];
    matchData.pt3D = vec_pts3D[cam];
    
    // create indMatch3D2D
    std::vector<pair<IndexT, IndexT> > indMatch3D2D;
    indMatch3D2D.reserve(matchData.pt2D.cols());
    const auto &occurrences = vec_occurrences[cam];
    for(const auto &ass : occurrences)
    {
      // recopy the associations IDs in the vector
      indMatch3D2D.push_back(ass.first);
    }
    
    vec_locResults.emplace_back(matchData, indMatch3D2D, pose, intrinsics, std::vector<voctree::DocMatch>(), refineOk);
  }
  
  if(!refineOk)
  {
    POPART_COUT("[poseEstimation]\tRefine failed.");
    return false;
  }
  
  return true;
  
}
  
#endif //HAVE_OPENGV

// subposes is n-1 as we consider the first camera as the main camera and the 
// reference frame of the grid
bool CCTagLocalizer::localizeRig_naive(const std::vector<std::unique_ptr<features::Regions> > & vec_queryRegions,
                                 const std::vector<std::pair<std::size_t, std::size_t> > &imageSize,
                                 const LocalizerParameters *parameters,
                                 std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                                 const std::vector<geometry::Pose3 > &vec_subPoses,
                                 geometry::Pose3 &rigPose,
                                 std::vector<LocalizationResult>& vec_localizationResults)
{
  const CCTagLocalizer::Parameters *param = static_cast<const CCTagLocalizer::Parameters *>(parameters);
  if(!param)
  {
    throw std::invalid_argument("The CCTag localizer parameters are not in the right format.");
  }

  const size_t numCams = vec_queryRegions.size();
  
  assert(numCams==vec_queryIntrinsics.size());
  assert(numCams==vec_subPoses.size()-1);
  assert(numCams==imageSize.size());

  vec_localizationResults.resize(numCams);
    
  // this is basic, just localize each camera alone
  //@todo parallelize?
  std::vector<bool> isLocalized(numCams, false);
  for(size_t i = 0; i < numCams; ++i)
  {
    isLocalized[i] = localize(vec_queryRegions[i], imageSize[i], param, true /*useInputIntrinsics*/, vec_queryIntrinsics[i], vec_localizationResults[i]);
    if(!isLocalized[i])
    {
      POPART_CERR("Could not localize camera " << i);
      // even if it is not localize we can try to go on and do with the cameras we have
    }
  }
  
  // ** 'easy' cases in which we don't need further processing **
  
  const std::size_t numLocalizedCam = std::count(isLocalized.begin(), isLocalized.end(), true);
  
  // no camera has be localized
  if(numLocalizedCam == 0)
  {
    POPART_COUT("No camera has been localized!!!");
    return false;
  }
  
  POPART_COUT("Localized cameras: " << numLocalizedCam << "/" << numCams);
  
  // if there is only one camera (the main one)
  if(numCams==1)
  {
    // there is only the main camera, not much else to do, the position is already
    // refined by the call to localize
    //set the pose
    rigPose = vec_localizationResults[0].getPose();
  }
  else
  {

    // find the index of the first localized camera
    const std::size_t idx = std::distance(isLocalized.begin(), 
                                          std::find(isLocalized.begin(), isLocalized.end(), true));
    
    // useless safeguard as there should be at least 1 element at this point but
    // better safe than sorry
    assert(idx < isLocalized.size());
    
    POPART_COUT("Index of the first localized camera: " << idx);
    // if the only localized camera is the main camera
    if(idx==0)
    {
      // just give its pose
      rigPose = vec_localizationResults[0].getPose();
    }
    else
    {
      // main camera: q1 ~ [R1 t1] Q = [I 0] A   where A = [R1 t1] Q  
      // another camera: q2 ~ [R2 t2] Q = [R2 t2]*inv([R1 t1]) A   and subPose12 = [R12 t12] = [R2 t2]*inv([R1 t1])
      // with the localization localize() we have computed [R2 t2], hence:
      // q2 ~ [R2 t2] Q = [R12 t12]*inv([R12 t12]) * [R2 t2] Q
      // and inv([R12 t12]) * [R2 t2] is the pose of the main camera
      
      // recover the rig pose using the subposes
      rigPose = vec_subPoses[idx-1].inverse() * vec_localizationResults[idx].getPose();
    }
  }
  
  // ** otherwise run a BA with the localized cameras
  const bool refineOk = refineRigPose(vec_subPoses, vec_localizationResults, rigPose);
  
  if(!refineOk)
  {
    POPART_COUT("[poseEstimation]\tRig pose refinement failed.");
    return false;
  }
  
  updateRigPoses(vec_localizationResults, rigPose, vec_subPoses);
  
  return true;
}


void CCTagLocalizer::getAllAssociations(const features::CCTAG_Regions &queryRegions,
                                        const std::pair<std::size_t, std::size_t> &imageSize,
                                        const CCTagLocalizer::Parameters &param,
                                        std::map< std::pair<IndexT, IndexT>, std::size_t > &occurences, 
                                        Mat &pt2D,
                                        Mat &pt3D,
                                        const std::string& imagePath) const
{
  std::vector<IndexT> nearestKeyFrames;
  nearestKeyFrames.reserve(param._nNearestKeyFrames);
  
  kNearestKeyFrames(queryRegions,
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
    POPART_COUT("matching]\tFound "<< vec_featureMatches.size() <<" matches.");
    
    if(!param._visualDebug.empty() && !imagePath.empty())
    {
      namespace bfs = boost::filesystem;
      const sfm::View *mview = _sfm_data.GetViews().at(indexKeyFrame).get();
      const std::string queryImage = bfs::path(imagePath).stem().string();
      const std::string matchedImage = bfs::path(mview->s_Img_path).stem().string();
      const std::string matchedPath = (bfs::path(_sfm_data.s_root_path) /  bfs::path(mview->s_Img_path)).string();

      // the directory where to save the feature matches
      const auto baseDir = bfs::path(param._visualDebug) / queryImage;
      if((!bfs::exists(baseDir)))
      {
        POPART_COUT("created " << baseDir.string());
        bfs::create_directories(baseDir);
      }
      
      // the final filename for the output svg file as a composition of the query
      // image and the matched image
      auto outputName = baseDir / queryImage;
      outputName += "_";
      outputName += matchedImage;
      outputName += ".svg";
      
      const bool showNotMatched = true;
      features::saveCCTagMatches2SVG(imagePath, 
                                     imageSize, 
                                     queryRegions,
                                     matchedPath,
                                     std::make_pair(mview->ui_width, mview->ui_height), 
                                     _regions_per_view.at(indexKeyFrame)._regions,
                                     vec_featureMatches,
                                     outputName.string(),
                                     showNotMatched ); 
    }
    
    // Recover the 2D-3D associations from the matches 
    // Each matched feature in the current similar image is associated to a 3D point,
    // hence we can recover the 2D-3D associations to estimate the pose
    
    // Get the 3D points associated to each matched feature
    for(const matching::IndMatch& featureMatch : vec_featureMatches)
    {
      // the ID of the 3D point
      const IndexT pt3D_id = matchedRegions._associated3dPoint[featureMatch._j];
      const IndexT pt2D_id = featureMatch._i;
      
      const auto key = std::make_pair(pt3D_id, pt2D_id);
      if(occurences.count(key))
      {
        occurences[key]++;
      }
      else
      {
        occurences[key] = 1;
      }
    }
  }
      
  const size_t numCollectedPts = occurences.size();
  POPART_COUT("[matching]\tCollected "<< numCollectedPts <<" associations.");
  
  {
    // just debugging statistics, this block can be safely removed    
    std::size_t maxOcc = 0;
    for(const auto &idx : occurences)
    {
      const auto &key = idx.first;
      const auto &value = idx.second;
       POPART_COUT("[matching]\tAssociations "
               << key.first << "," << key.second <<"] found " 
               << value << " times.");
       if(value > maxOcc)
         maxOcc = value;
    }
    
    std::size_t numOccTreated = 0;
    for(std::size_t value = 1; value < maxOcc; ++value)
    {
      std::size_t counter = 0;
      for(const auto &idx : occurences)
      {
        if(idx.second == value)
        {
          ++counter;
        }
      }
      if(counter>0)
        POPART_COUT("[matching]\tThere are " << counter
                    << " associations occurred " << value << " times ("
                    << 100.0 * counter / (double) numCollectedPts << "%)");
      numOccTreated += counter;
      if(numOccTreated >= numCollectedPts)
        break;
    }
  }

  pt2D = Mat2X(2, numCollectedPts);
  pt3D = Mat3X(3, numCollectedPts);
      
  size_t index = 0;
  for(const auto &idx : occurences)
  {
    // recopy all the points in the matching structure
    const IndexT pt3D_id = idx.first.first;
    const IndexT pt2D_id = idx.first.second;
      
    pt2D.col(index) = queryRegions.GetRegionPosition(pt2D_id);
    pt3D.col(index) = _sfm_data.GetLandmarks().at(pt3D_id).X;
      ++index;
  }
}





void kNearestKeyFrames(const features::CCTAG_Regions & queryRegions,
                       const CCTagRegionsPerViews & regionsPerView,
                       std::size_t nNearestKeyFrames,
                       std::vector<IndexT> & kNearestFrames,
                       const float similarityThreshold /*=.0f*/)
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
  kNearestFrames.reserve(nNearestKeyFrames);
  for (auto rit = sortedViewSimilarities.crbegin(); rit != sortedViewSimilarities.crend(); ++rit)
  {
    if(rit->first < similarityThreshold)
      // since it is ordered, the first having smaller similarity guarantees that
      // there won't be other useful kframes
      break;
    
    kNearestFrames.push_back(rit->second);
    ++counter;
    
    if (counter == nNearestKeyFrames)
      break;
  }
}
 
void viewMatching(const features::CCTAG_Regions & regionsA,
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
 
 
 
float viewSimilarity(const features::CCTAG_Regions & regionsA,
                     const features::CCTAG_Regions & regionsB)
{
  assert(regionsA.DescriptorLength() == regionsB.DescriptorLength()); 
  
  const std::bitset<128> descriptorViewA = constructCCTagViewDescriptor(regionsA.Descriptors());
  const std::bitset<128> descriptorViewB = constructCCTagViewDescriptor(regionsB.Descriptors());
  
  // The similarity is the sum of all the cctags sharing the same id visible in both views.
  return (descriptorViewA & descriptorViewB).count();
}

std::bitset<128> constructCCTagViewDescriptor(const std::vector<CCTagDescriptor> & vCCTagDescriptors)
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

