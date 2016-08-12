#include "VoctreeLocalizer.hpp"
#include "rigResection.hpp"
#include "optimization.hpp"

#include <openMVG/sfm/sfm_data_io.hpp>
#include <openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp>
#include <openMVG/sfm/sfm_data_BA_ceres.hpp>
#include <openMVG/features/io_regions_type.hpp>
#include <openMVG/features/svgVisualization.hpp>
#ifdef HAVE_CCTAG
#include <openMVG/features/cctag/SIFT_CCTAG_describer.hpp>
#endif
#include <nonFree/sift/SIFT_float_describer.hpp>
#include <openMVG/matching/regions_matcher.hpp>
#include <openMVG/matching_image_collection/Matcher.hpp>
#include <openMVG/matching/matcher_kdtree_flann.hpp>
#include <openMVG/matching_image_collection/F_ACRobust.hpp>
#include <openMVG/numeric/numeric.h>
#include <openMVG/robust_estimation/guided_matching.hpp>
#include <openMVG/logger.hpp>

#include <third_party/progress/progress.hpp>
//#include <cereal/archives/json.hpp>

#include <boost/filesystem.hpp>

#include <algorithm>
#include <chrono>

namespace openMVG {
namespace localization {

std::ostream& operator<<( std::ostream& os, const voctree::Document &doc )	
{
  os << "[ ";
  for( const voctree::Word &w : doc )
  {
          os << w << ", ";
  }
  os << "];\n";
  return os;
}

std::ostream& operator<<(std::ostream& os, VoctreeLocalizer::Algorithm a)
{
  switch(a)
  {
  case VoctreeLocalizer::Algorithm::FirstBest: os << "FirstBest";
    break;
  case VoctreeLocalizer::Algorithm::BestResult: os << "BestResult";
    break;
  case VoctreeLocalizer::Algorithm::AllResults: os << "AllResults";
    break;
  case VoctreeLocalizer::Algorithm::Cluster: os << "Cluster";
    break;
  default: 
    os << "Unknown algorithm!";
    throw std::invalid_argument("Unrecognized algorithm!");
  }
  return os;
}

std::istream& operator>>(std::istream &in, VoctreeLocalizer::Algorithm &a)
{
  POPART_COUT("operator>>");
  int i;
  in >> i;
  a = static_cast<VoctreeLocalizer::Algorithm> (i);
  return in;
}	

VoctreeLocalizer::Algorithm VoctreeLocalizer::initFromString(const std::string &value)
{
  if(value=="FirstBest")
    return VoctreeLocalizer::Algorithm::FirstBest;
  else if(value=="AllResults")
    return VoctreeLocalizer::Algorithm::AllResults;
  else if(value=="BestResult")
    throw std::invalid_argument("BestResult not yet implemented");
  else if(value=="Cluster")
    throw std::invalid_argument("Cluster not yet implemented");
  else
    throw std::invalid_argument("Unrecognized algorithm \"" + value + "\"!");
}


VoctreeLocalizer::VoctreeLocalizer(const std::string &sfmFilePath,
                                   const std::string &descriptorsFolder,
                                   const std::string &vocTreeFilepath,
                                   const std::string &weightsFilepath
#ifdef HAVE_CCTAG
                                   , bool useSIFT_CCTAG
#endif
                                  ) : ILocalizer()
{
  using namespace openMVG::features;
  // init the feature extractor
#ifdef HAVE_CCTAG
  if(useSIFT_CCTAG)
  {
    POPART_COUT("SIFT_CCTAG_Image_describer");
    _image_describer = new features::SIFT_CCTAG_Image_describer();  
  }
  else
  {
#if USE_SIFT_FLOAT
    POPART_COUT("SIFT_float_describer");
    _image_describer = new features::SIFT_float_describer();
#else
    POPART_COUT("SIFT_Image_describer");
    _image_describer = new features::SIFT_Image_describer();
#endif
  }
#else //HAVE_CCTAG
#if USE_SIFT_FLOAT
    POPART_COUT("SIFT_float_describer");
    _image_describer = new features::SIFT_float_describer();
#else
    POPART_COUT("SIFT_Image_describer");
    _image_describer = new features::SIFT_Image_describer();
#endif
#endif //HAVE_CCTAG
  
  // load the sfm data containing the 3D reconstruction info
  POPART_COUT("Loading SFM data...");
  if (!Load(_sfm_data, sfmFilePath, sfm::ESfM_Data::ALL)) 
  {
    POPART_CERR("The input SfM_Data file " << sfmFilePath << " cannot be read!");
    POPART_CERR("\n\nIf the error says \"JSON Parsing failed - provided NVP not found\" "
            "it's likely that you have to convert your sfm_data to a recent version supporting "
            "polymorphic Views. You can run the python script convertSfmData.py to update an existing sfmdata.");
    _isInit = false;
    return;
  }

  POPART_COUT("SfM data loaded from " << sfmFilePath << " containing: ");
  POPART_COUT("\tnumber of views      : " << _sfm_data.GetViews().size());
  POPART_COUT("\tnumber of poses      : " << _sfm_data.GetPoses().size());
  POPART_COUT("\tnumber of points     : " << _sfm_data.GetLandmarks().size());
  POPART_COUT("\tnumber of intrinsics : " << _sfm_data.GetIntrinsics().size());

  // load the features and descriptors
  // initially we need all the feature in order to create the database
  // then we can store only those associated to 3D points
  //? can we use Feature_Provider to load the features and filter them later?

  _isInit = initDatabase(vocTreeFilepath, weightsFilepath, descriptorsFolder);
}

bool VoctreeLocalizer::localize(const std::unique_ptr<features::Regions> &genQueryRegions,
                                const std::pair<std::size_t, std::size_t> &imageSize,
                                const LocalizerParameters *param,
                                bool useInputIntrinsics,
                                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                                LocalizationResult & localizationResult,
                                const std::string& imagePath)
{
  const Parameters *voctreeParam = static_cast<const Parameters *>(param);
  if(!voctreeParam)
  {
    // error!
    throw std::invalid_argument("The parameters are not in the right format!!");
  }

#if USE_SIFT_FLOAT  
  features::SIFT_Float_Regions *queryRegions = dynamic_cast<features::SIFT_Float_Regions*> (genQueryRegions.get());
#else
  features::SIFT_Regions *queryRegions = dynamic_cast<features::SIFT_Regions*> (genQueryRegions.get());
#endif
  
  if(!queryRegions)
  {
    // error, the casting did not work
    throw std::invalid_argument("Error while converting the input Regions. Only SIFT regions are allowed for the voctree localizer.");
  }
  
  switch(voctreeParam->_algorithm)
  {
    case Algorithm::FirstBest:
    return localizeFirstBestResult(*queryRegions,
                                   imageSize,
                                   *voctreeParam,
                                   useInputIntrinsics,
                                   queryIntrinsics,
                                   localizationResult,
                                   imagePath);
    case Algorithm::BestResult: throw std::invalid_argument("BestResult not yet implemented");
    case Algorithm::AllResults:
    return localizeAllResults(*queryRegions,
                              imageSize,
                              *voctreeParam,
                              useInputIntrinsics,
                              queryIntrinsics,
                              localizationResult,
                              imagePath);
    case Algorithm::Cluster: throw std::invalid_argument("Cluster not yet implemented");
    default: throw std::invalid_argument("Unknown algorithm type");
  }
}

bool VoctreeLocalizer::localize(const image::Image<unsigned char> & imageGrey,
                                const LocalizerParameters *param,
                                bool useInputIntrinsics,
                                cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                                LocalizationResult &localizationResult,
                                const std::string& imagePath /* = std::string() */)
{
  // A. extract descriptors and features from image
  POPART_COUT("[features]\tExtract SIFT from query image");
#if USE_SIFT_FLOAT
  std::unique_ptr<features::Regions> tmpQueryRegions(new features::SIFT_Float_Regions());
#else
  std::unique_ptr<features::Regions> tmpQueryRegions(new features::SIFT_Regions());
#endif
  auto detect_start = std::chrono::steady_clock::now();
  _image_describer->Set_configuration_preset(param->_featurePreset);
  _image_describer->Describe(imageGrey, tmpQueryRegions, nullptr);
  auto detect_end = std::chrono::steady_clock::now();
  auto detect_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(detect_end - detect_start);
  POPART_COUT("[features]\tExtract SIFT done: found " << tmpQueryRegions->RegionCount() << " features in " << detect_elapsed.count() << " [ms]");

  const std::pair<std::size_t, std::size_t> queryImageSize = std::make_pair(imageGrey.Width(), imageGrey.Height());

  // if debugging is enable save the svg image with the extracted features
  if(!param->_visualDebug.empty() && !imagePath.empty())
  {
    namespace bfs = boost::filesystem;
    saveFeatures2SVG(imagePath,
                     queryImageSize,
                     tmpQueryRegions->GetRegionsPositions(),
                     param->_visualDebug + "/" + bfs::path(imagePath).stem().string() + ".svg");
  }

  return localize(tmpQueryRegions,
                  queryImageSize,
                  param,
                  useInputIntrinsics,
                  queryIntrinsics,
                  localizationResult,
                  imagePath);
}

//@fixme deprecated.. now inside initDatabase

bool VoctreeLocalizer::loadReconstructionDescriptors(const sfm::SfM_Data & sfm_data,
                                                     const std::string & feat_directory)
{
  C_Progress_display my_progress_bar(sfm_data.GetViews().size(),
                                     std::cout, "\n- Regions Loading -\n");

  POPART_COUT("Build observations per view");
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

  POPART_COUT("Load Features and Descriptors per view");
  // Read for each view the corresponding regions and store them
  for(sfm::Views::const_iterator iter = sfm_data.GetViews().begin();
          iter != sfm_data.GetViews().end(); ++iter, ++my_progress_bar)
  {
    const IndexT id_view = iter->second->id_view;
    Reconstructed_RegionsT& reconstructedRegion = _regions_per_view[id_view];

    const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
    const std::string basename = stlplus::basename_part(sImageName);
    const std::string featFilepath = stlplus::create_filespec(feat_directory, basename, ".feat");
    const std::string descFilepath = stlplus::create_filespec(feat_directory, basename, ".desc");
    //    std::cout << "Feat: " << featFilepath << std::endl;
    //    std::cout << "Desc: " << descFilepath << std::endl;

    if(!reconstructedRegion._regions.Load(featFilepath, descFilepath))
    {
      POPART_CERR("Invalid regions files for the view: " << sImageName);
      return false;
    }

    // Filter descriptors to keep only the 3D reconstructed points
    reconstructedRegion.filterRegions(observationsPerView[id_view]);
  }
  return true;
}

/**
 * @brief Initialize the database...
 */
bool VoctreeLocalizer::initDatabase(const std::string & vocTreeFilepath,
                                    const std::string & weightsFilepath,
                                    const std::string & feat_directory)
{

  bool withWeights = !weightsFilepath.empty();

  // Load vocabulary tree
  POPART_COUT("Loading vocabulary tree...");

  _voctree.load(vocTreeFilepath);
  POPART_COUT("tree loaded with" << endl << "\t" << _voctree.levels() << " levels" 
          << endl << "\t" << _voctree.splits() << " branching factor");

  POPART_COUT("Creating the database...");
  // Add each object (document) to the database
  _database = voctree::Database(_voctree.words());
  if(withWeights)
  {
    POPART_COUT("Loading weights...");
    _database.loadWeights(weightsFilepath);
  }
  else
  {
    POPART_COUT("No weights specified, skipping...");
  }
  
  // Load the descriptors and the features related to the images
  // for every image, pass the descriptors through the vocabulary tree and
  // add its visual words to the database.
  // then only store the feature and descriptors that have a 3D point associated
  POPART_COUT("Build observations per view");
  C_Progress_display my_progress_bar(_sfm_data.GetViews().size(),
                                     std::cout, "\n- Load Features and Descriptors per view -\n");

  // Build observations per view
  std::map<IndexT, std::vector<FeatureInImage> > observationsPerView;
  for(auto landmarkValue : _sfm_data.structure)
  {
    IndexT trackId = landmarkValue.first;
    sfm::Landmark& landmark = landmarkValue.second;
    for(auto obs : landmark.obs)
    {
      const IndexT viewId = obs.first;
      const sfm::Observation& obs2d = obs.second;
      observationsPerView[viewId].emplace_back(obs2d.id_feat, trackId);
    }
  }
  for(auto featuresInImage : observationsPerView)
  {
    std::sort(featuresInImage.second.begin(), featuresInImage.second.end());
  }

  // Read for each view the corresponding Regions and store them
  for(const auto &iter : _sfm_data.GetViews())
  {
    const std::shared_ptr<sfm::View> currView = iter.second;
    const IndexT id_view = currView->id_view;
    Reconstructed_RegionsT& currRecoRegions = _regions_per_view[id_view];

    const std::string sImageName = stlplus::create_filespec(_sfm_data.s_root_path, currView.get()->s_Img_path);
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

    if(!currRecoRegions._regions.Load(featFilepath, descFilepath))
    {
      std::cerr << "Invalid regions files for the view: " << sImageName << std::endl;
      return false;
    }
    
    voctree::SparseHistogram histo;
    std::vector<voctree::Word> words = _voctree.quantize(currRecoRegions._regions.Descriptors());
    voctree::computeSparseHistogram(words, histo);
    _database.insert(id_view, histo);

    // Filter descriptors to keep only the 3D reconstructed points
    currRecoRegions.filterRegions(observationsPerView[id_view]);
    ++my_progress_bar;
  }
  
  return true;
}

bool VoctreeLocalizer::localizeFirstBestResult(const features::SIFT_Regions &queryRegions,
                                               const std::pair<std::size_t, std::size_t> queryImageSize,
                                               const Parameters &param,
                                               bool useInputIntrinsics,
                                               cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                                               LocalizationResult &localizationResult,
                                               const std::string& imagePath /*= std::string()*/)
{
  // A. Find the (visually) similar images in the database 
  POPART_COUT("[database]\tRequest closest images from voctree");
  // pass the descriptors through the vocabulary tree to get the visual words
  // associated to each feature
  std::vector<voctree::Word> requestImageWords = _voctree.quantize(queryRegions.Descriptors());
  
  // Request closest images from voctree
  std::vector<voctree::DocMatch> matchedImages;
  _database.find(requestImageWords, param._numResults, matchedImages);
  
//  // just debugging bla bla
//  // for each similar image found print score and number of features
//  for(const voctree::DocMatch & currMatch : matchedImages )
//  {
//    // get the corresponding index of the view
//    const IndexT matchedViewIndex = currMatch.id;
//    // get the view handle
//    const std::shared_ptr<sfm::View> matchedView = _sfm_data.views[matchedViewIndex];
//    POPART_COUT( "[database]\t\t match " << matchedView->s_Img_path 
//            << " [docid: "<< currMatch.id << "]"
//            << " with score " << currMatch.score 
//            << " and it has "  << _regions_per_view[matchedViewIndex]._regions.RegionCount() 
//            << " features with 3D points");
//  }

  POPART_COUT("[matching]\tBuilding the matcher");
  matching::RegionsMatcherT<MatcherT> matcher(queryRegions);
  
  sfm::Image_Localizer_Match_Data resectionData;
  std::vector<pair<IndexT, IndexT> > associationIDs;
  geometry::Pose3 pose;
 
  // B. for each found similar image, try to find the correspondences between the 
  // query image and the similar image
  for(const voctree::DocMatch& matchedImage : matchedImages)
  {
    // minimum number of points that allows a reliable 3D reconstruction
    const size_t minNum3DPoints = 5;
    
    // the view index of the current matched image
    const IndexT matchedViewIndex = matchedImage.id;
    // the handler to the current view
    const std::shared_ptr<sfm::View> matchedView = _sfm_data.views[matchedViewIndex];
    
    // safeguard: we should match the query image with an image that has at least
    // some 3D points visible --> if it has 0 3d points it is likely that it is an
    // image of the dataset that was not reconstructed
    if(_regions_per_view[matchedViewIndex]._regions.RegionCount() < minNum3DPoints)
    {
      POPART_COUT("[matching]\tSkipping matching with " << matchedView->s_Img_path << " as it has too few visible 3D points (" << _regions_per_view[matchedViewIndex]._regions.RegionCount() << ")");
      continue;
    }
    else
    {
      POPART_COUT("[matching]\tTrying to match the query image with " << matchedView->s_Img_path);
    }
    
    // its associated reconstructed regions
    const Reconstructed_RegionsT& matchedRegions = _regions_per_view[matchedViewIndex];
    // its associated intrinsics
    // this is just ugly!
    const cameras::IntrinsicBase *matchedIntrinsicsBase = _sfm_data.intrinsics[matchedView->id_intrinsic].get();
    if ( !isPinhole(matchedIntrinsicsBase->getType()) )
    {
      //@fixme maybe better to throw something here
      POPART_CERR("Only Pinhole cameras are supported!");
      return false;
    }
    const cameras::Pinhole_Intrinsic *matchedIntrinsics = (const cameras::Pinhole_Intrinsic*)(matchedIntrinsicsBase);
     
    std::vector<matching::IndMatch> vec_featureMatches;
    bool matchWorked = robustMatching( matcher, 
                                      // pass the input intrinsic if they are valid, null otherwise
                                      (useInputIntrinsics) ? &queryIntrinsics : nullptr,
                                      matchedRegions,
                                      matchedIntrinsics,
                                      param._fDistRatio,
                                      param._matchingError,
                                      param._useGuidedMatching,
                                      queryImageSize,
                                      std::make_pair(matchedView->ui_width, matchedView->ui_height), 
                                      vec_featureMatches,
                                      param._matchingEstimator);
    if (!matchWorked)
    {
      POPART_COUT("[matching]\tMatching with " << matchedView->s_Img_path << " failed! Skipping image");
      continue;
    }
    else
    {
      POPART_COUT("[matching]\tFound " << vec_featureMatches.size() << " matches");
    }
    assert(vec_featureMatches.size()>0);
    
    if(!param._visualDebug.empty() && !imagePath.empty())
    {
      namespace bfs = boost::filesystem;
      const sfm::View *mview = _sfm_data.GetViews().at(matchedViewIndex).get();
      const std::string queryimage = bfs::path(imagePath).stem().string();
      const std::string matchedImage = bfs::path(mview->s_Img_path).stem().string();
      const std::string matchedPath = (bfs::path(_sfm_data.s_root_path) /  bfs::path(mview->s_Img_path)).string();
      
      features::saveMatches2SVG(imagePath,
                      queryImageSize,
                      queryRegions.GetRegionsPositions(),
                      matchedPath,
                      std::make_pair(mview->ui_width, mview->ui_height),
                      _regions_per_view[matchedViewIndex]._regions.GetRegionsPositions(),
                      vec_featureMatches,
                      param._visualDebug + "/" + queryimage + "_" + matchedImage + ".svg"); 
    }
    
    // C. recover the 2D-3D associations from the matches 
    // Each matched feature in the current similar image is associated to a 3D point,
    // hence we can recover the 2D-3D associations to estimate the pose
    // Prepare data for resection
    resectionData = sfm::Image_Localizer_Match_Data();
    resectionData.pt2D = Mat2X(2, vec_featureMatches.size());
    resectionData.pt3D = Mat3X(3, vec_featureMatches.size());
    associationIDs.clear();
    associationIDs.reserve(vec_featureMatches.size());

    // Get the 3D points associated to each matched feature
    std::size_t index = 0;
    for(const matching::IndMatch& featureMatch : vec_featureMatches)
    {
      assert(vec_featureMatches.size()>index);
      // the ID of the 3D point
      const IndexT trackId3D = matchedRegions._associated3dPoint[featureMatch._j];

      // prepare data for resectioning
      resectionData.pt3D.col(index) = _sfm_data.GetLandmarks().at(trackId3D).X;

      const Vec2 feat = queryRegions.GetRegionPosition(featureMatch._i);
      resectionData.pt2D.col(index) = feat;
      
      associationIDs.emplace_back(trackId3D, featureMatch._i);

      ++index;
    }
    // estimate the pose
    // Do the resectioning: compute the camera pose.
    resectionData.error_max = param._errorMax;
    POPART_COUT("[poseEstimation]\tEstimating camera pose...");
    bool bResection = sfm::SfM_Localizer::Localize(queryImageSize,
                                                   // pass the input intrinsic if they are valid, null otherwise
                                                   (useInputIntrinsics) ? &queryIntrinsics : nullptr,
                                                   resectionData,
                                                   pose,
                                                   param._resectionEstimator);

    if(!bResection)
    {
      POPART_COUT("[poseEstimation]\tResection FAILED");
      continue;
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
      POPART_COUT("K estimated\n" << K_);
      queryIntrinsics.setK(K_);
      queryIntrinsics.setWidth(queryImageSize.first);
      queryIntrinsics.setHeight(queryImageSize.second);
    }

    // D. refine the estimated pose
    POPART_COUT("[poseEstimation]\tRefining estimated pose");
    bool refineStatus = sfm::SfM_Localizer::RefinePose(&queryIntrinsics, 
                                                       pose, 
                                                       resectionData, 
                                                       true /*b_refine_pose*/, 
                                                       param._refineIntrinsics /*b_refine_intrinsic*/);
    if(!refineStatus)
      POPART_COUT("[poseEstimation]\tRefine pose could not improve the estimation of the camera pose.");
    
    {
      // just temporary code to evaluate the estimated pose @todo remove it
      const geometry::Pose3 &referencePose = _sfm_data.poses[matchedViewIndex];
      POPART_COUT("R refined\n" << pose.rotation());
      POPART_COUT("t refined\n" << pose.translation());
      POPART_COUT("K refined\n" << queryIntrinsics.K());
      POPART_COUT("R_gt\n" << referencePose.rotation());
      POPART_COUT("t_gt\n" << referencePose.translation());
      POPART_COUT("angular difference: " << R2D(getRotationMagnitude(pose.rotation()*referencePose.rotation().inverse())) << "deg");
      POPART_COUT("center difference: " << (pose.center()-referencePose.center()).norm());
      POPART_COUT("err = [err; " << R2D(getRotationMagnitude(pose.rotation()*referencePose.rotation().inverse())) << ", "<< (pose.center()-referencePose.center()).norm() << "];");
    }
    localizationResult = LocalizationResult(resectionData, associationIDs, pose, queryIntrinsics, matchedImages, true);
    break;
  }
  //@todo deal with unsuccesful case...
  return localizationResult.isValid();
  
 } 

bool VoctreeLocalizer::localizeAllResults(const features::SIFT_Regions &queryRegions,
                                          const std::pair<std::size_t, std::size_t> queryImageSize,
                                          const Parameters &param,
                                          bool useInputIntrinsics,
                                          cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                                          LocalizationResult &localizationResult,
                                          const std::string& imagePath)
{
  
  sfm::Image_Localizer_Match_Data resectionData;
  std::map< std::pair<IndexT, IndexT>, std::size_t > occurences;
  
  // get all the association from the database images
  std::vector<voctree::DocMatch> matchedImages;
  getAllAssociations(queryRegions,
                     queryImageSize,
                     param,
                     useInputIntrinsics,
                     queryIntrinsics,
                     occurences,
                     resectionData.pt2D,
                     resectionData.pt3D,
                     matchedImages,
                     imagePath);

  const std::size_t numCollectedPts = occurences.size();
  std::vector<pair<IndexT, IndexT> > associationIDs;
  associationIDs.reserve(numCollectedPts);

  for(const auto &ass : occurences)
  {
    // recopy the associations IDs in the vector
    associationIDs.push_back(ass.first);
  }

  geometry::Pose3 pose;
  
  // estimate the pose
  // Do the resectioning: compute the camera pose.
  resectionData.error_max = param._errorMax;
  POPART_COUT("[poseEstimation]\tEstimating camera pose...");
  const bool bResection = sfm::SfM_Localizer::Localize(queryImageSize,
                                                      // pass the input intrinsic if they are valid, null otherwise
                                                      (useInputIntrinsics) ? &queryIntrinsics : nullptr,
                                                      resectionData,
                                                      pose,
                                                      param._resectionEstimator);

  if(!bResection)
  {
    POPART_COUT("[poseEstimation]\tResection FAILED");
    if(!param._visualDebug.empty() && !imagePath.empty())
    {
      namespace bfs = boost::filesystem;
      features::saveFeatures2SVG(imagePath, 
                       queryImageSize, 
                         resectionData.pt2D,
                         param._visualDebug + "/" + bfs::path(imagePath).stem().string() + ".associations.svg");
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
    queryIntrinsics.setWidth(queryImageSize.first);
    queryIntrinsics.setHeight(queryImageSize.second);
  }

  // E. refine the estimated pose
  POPART_COUT("[poseEstimation]\tRefining estimated pose");
  bool refineStatus = sfm::SfM_Localizer::RefinePose(&queryIntrinsics,
                                                     pose,
                                                     resectionData,
                                                     true /*b_refine_pose*/,
                                                     param._refineIntrinsics /*b_refine_intrinsic*/);
  if(!refineStatus)
    POPART_COUT("Refine pose could not improve the estimation of the camera pose.");

  if(!param._visualDebug.empty() && !imagePath.empty())
  {
    namespace bfs = boost::filesystem;
    features::saveFeatures2SVG(imagePath,
                     queryImageSize,
                     resectionData.pt2D,
                     param._visualDebug + "/" + bfs::path(imagePath).stem().string() + ".associations.svg",
                     &resectionData.vec_inliers);
  }

  localizationResult = LocalizationResult(resectionData, associationIDs, pose, queryIntrinsics, matchedImages, true);

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

void VoctreeLocalizer::getAllAssociations(const features::SIFT_Regions &queryRegions,
                                          const std::pair<std::size_t, std::size_t> &imageSize,
                                          const Parameters &param,
                                          bool useInputIntrinsics,
                                          const cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics,
                                          std::map< std::pair<IndexT, IndexT>, std::size_t > &occurences,
                                          Mat &pt2D,
                                          Mat &pt3D,
                                          std::vector<voctree::DocMatch>& matchedImages,
                                          const std::string& imagePath) const
{
  // A. Find the (visually) similar images in the database 
  // pass the descriptors through the vocabulary tree to get the visual words
  // associated to each feature
  POPART_COUT("[database]\tRequest closest images from voctree");
  std::vector<voctree::Word> requestImageWords = _voctree.quantize(queryRegions.Descriptors());
  
  // Request closest images from voctree
  _database.find(requestImageWords, (param._numResults==0) ? (_database.size()) : (param._numResults) , matchedImages);

//  // just debugging bla bla
//  // for each similar image found print score and number of features
//  for(const voctree::DocMatch& currMatch : matchedImages )
//  {
//    // get the view handle
//    const std::shared_ptr<sfm::View> matchedView = _sfm_data.views[currMatch.id];
//    POPART_COUT( "[database]\t\t match " << matchedView->s_Img_path 
//            << " [docid: "<< currMatch.id << "]"
//            << " with score " << currMatch.score 
//            << " and it has "  << _regions_per_view[currMatch.id]._regions.RegionCount() 
//            << " features with 3D points");
//  }


  POPART_COUT("[matching]\tBuilding the matcher");
  matching::RegionsMatcherT<MatcherT> matcher(queryRegions);

  
  std::map< std::pair<IndexT, IndexT>, std::size_t > repeated;
  
  // B. for each found similar image, try to find the correspondences between the 
  // query image adn the similar image
  // stop when param._maxResults successful matches have been found
  std::size_t goodMatches = 0;
  for(const voctree::DocMatch& matchedImage : matchedImages)
  {
    // minimum number of points that allows a reliable 3D reconstruction
    const size_t minNum3DPoints = 5;
    
    // the handler to the current view
    const std::shared_ptr<sfm::View> matchedView = _sfm_data.views.at(matchedImage.id);
    // its associated reconstructed regions
    const Reconstructed_RegionsT& matchedRegions = _regions_per_view.at(matchedImage.id);
    
    // safeguard: we should match the query image with an image that has at least
    // some 3D points visible --> if this is not true it is likely that it is an
    // image of the dataset that was not reconstructed
    if(matchedRegions._regions.RegionCount() < minNum3DPoints)
    {
      POPART_COUT("[matching]\tSkipping matching with " << matchedView->s_Img_path << " as it has too few visible 3D points");
      continue;
    }
    POPART_COUT("[matching]\tTrying to match the query image with " << matchedView->s_Img_path);
    
    // its associated intrinsics
    // this is just ugly!
    const cameras::IntrinsicBase *matchedIntrinsicsBase = _sfm_data.intrinsics.at(matchedView->id_intrinsic).get();
    if ( !isPinhole(matchedIntrinsicsBase->getType()) )
    {
      //@fixme maybe better to throw something here
      POPART_CERR("Only Pinhole cameras are supported!");
      return;
    }
    const cameras::Pinhole_Intrinsic *matchedIntrinsics = (const cameras::Pinhole_Intrinsic*)(matchedIntrinsicsBase);
     
    std::vector<matching::IndMatch> vec_featureMatches;
    bool matchWorked = robustMatching( matcher, 
                                      // pass the input intrinsic if they are valid, null otherwise
                                      (useInputIntrinsics) ? &queryIntrinsics : nullptr,
                                      matchedRegions,
                                      matchedIntrinsics,
                                      param._fDistRatio,
                                      param._matchingError,
                                      param._useGuidedMatching,
                                      imageSize,
                                      std::make_pair(matchedView->ui_width, matchedView->ui_height), 
                                      vec_featureMatches,
                                      param._matchingEstimator);
    if (!matchWorked)
    {
//      POPART_COUT("[matching]\tMatching with " << matchedView->s_Img_path << " failed! Skipping image");
      continue;
    }
    else
    {
      POPART_COUT("[matching]\tFound " << vec_featureMatches.size() << " matches");
    }
    assert(vec_featureMatches.size()>0);
    
    // if debug is enable save the matches between the query image and the current matching image
    // It saves the feature matches in a folder with the same name as the query
    // image, if it does not exist it will create it. The final svg file will have
    // a name like this: queryImage_matchedImage.svg placed in the following directory:
    // param._visualDebug/queryImage/
    if(!param._visualDebug.empty() && !imagePath.empty())
    {
      namespace bfs = boost::filesystem;
      const auto &matchedViewIndex = matchedImage.id;
      const sfm::View *mview = _sfm_data.GetViews().at(matchedViewIndex).get();
      // the current query image without extension
      const auto queryImage = bfs::path(imagePath).stem();
      // the matching image without extension
      const auto matchedImage = bfs::path(mview->s_Img_path).stem();
      // the full path of the matching image
      const auto matchedPath = (bfs::path(_sfm_data.s_root_path) /  bfs::path(mview->s_Img_path)).string();

      // the directory where to save the feature matches
      const auto baseDir = bfs::path(param._visualDebug) / queryImage;
      if((!bfs::exists(baseDir)))
      {
        POPART_COUT("created " << baseDir.string());
        bfs::create_directories(baseDir);
      }
      
      // damn you, boost, what does it take to make the operator "+"?
      // the final filename for the output svg file as a composition of the query
      // image and the matched image
      auto outputName = baseDir / queryImage;
      outputName += "_";
      outputName += matchedImage;
      outputName += ".svg";

      features::saveMatches2SVG(imagePath,
                                imageSize,
                                queryRegions.GetRegionsPositions(),
                                matchedPath,
                                std::make_pair(mview->ui_width, mview->ui_height),
                                _regions_per_view.at(matchedViewIndex)._regions.GetRegionsPositions(),
                                vec_featureMatches,
                                outputName.string()); 
    }
    
    // C. recover the 2D-3D associations from the matches 
    // Each matched feature in the current similar image is associated to a 3D point
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
    ++goodMatches;
    if((param._maxResults !=0) && (goodMatches == param._maxResults))
    { 
      // let's say we have enough features
      POPART_COUT("[matching]\tgot enough point from " << param._maxResults << " images");
      break;
    }
  }
  
  const size_t numCollectedPts = occurences.size();
  
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

bool VoctreeLocalizer::robustMatching(matching::RegionsMatcherT<MatcherT> & matcher, 
                                      const cameras::IntrinsicBase * queryIntrinsicsBase,   // the intrinsics of the image we are using as reference
                                      const Reconstructed_RegionsT & matchedRegions,
                                      const cameras::IntrinsicBase * matchedIntrinsicsBase,
                                      const float fDistRatio,
                                      const double matchingError,
                                      const bool b_guided_matching,
                                      const std::pair<size_t,size_t> & imageSizeI,     // size of the first image @fixme change the API of the kernel!! 
                                      const std::pair<size_t,size_t> & imageSizeJ,     // size of the first image
                                      std::vector<matching::IndMatch> & vec_featureMatches,
                                      robust::EROBUST_ESTIMATOR estimator /*= robust::ROBUST_ESTIMATOR_ACRANSAC*/) const
{
  // get the intrinsics of the query camera
  if ((queryIntrinsicsBase != nullptr) && !isPinhole(queryIntrinsicsBase->getType()))
  {
    //@fixme maybe better to throw something here
    POPART_CERR("[matching]\tOnly Pinhole cameras are supported!");
    return false;
  }
  const cameras::Pinhole_Intrinsic *queryIntrinsics = (const cameras::Pinhole_Intrinsic*)(queryIntrinsicsBase);
  
  // get the intrinsics of the matched view
  if ((matchedIntrinsicsBase != nullptr) &&  !isPinhole(matchedIntrinsicsBase->getType()) )
  {
    //@fixme maybe better to throw something here
    POPART_CERR("[matching]\tOnly Pinhole cameras are supported!");
    return false;
  }
  const cameras::Pinhole_Intrinsic *matchedIntrinsics = (const cameras::Pinhole_Intrinsic*)(matchedIntrinsicsBase);
  
  const bool canBeUndistorted = (queryIntrinsicsBase != nullptr) && (matchedIntrinsicsBase != nullptr);
    
  // A. Putative Features Matching
  const bool matchWorked = matcher.Match(fDistRatio, matchedRegions._regions, vec_featureMatches);
  if (!matchWorked)
  {
    POPART_COUT("[matching]\tPutative matching failed.");
    return false;
  }
  assert(vec_featureMatches.size()>0);
  // prepare the data for geometric filtering: for each matched pair of features,
  // store them in two matrices
  Mat featuresI(2, vec_featureMatches.size());
  Mat featuresJ(2, vec_featureMatches.size());
  // fill the matrices with the features according to vec_featureMatches
  for(int i = 0; i < vec_featureMatches.size(); ++i)
  {
    const matching::IndMatch& match = vec_featureMatches[i];
    const Vec2 &queryPoint = matcher.getDatabaseRegions()->GetRegionPosition(match._i);
    const Vec2 &matchedPoint = matchedRegions._regions.GetRegionPosition(match._j);
    
    if(canBeUndistorted)
    {
      // undistort the points for the query image
      featuresI.col(i) = queryIntrinsics->get_ud_pixel(queryPoint);
      // undistort the points for the query image
      featuresJ.col(i) = matchedIntrinsics->get_ud_pixel(matchedPoint);
    }
    else
    {
      featuresI.col(i) = queryPoint;
      featuresJ.col(i) = matchedPoint;
    }
  }
  // perform the geometric filtering
  matching_image_collection::GeometricFilter_FMatrix_AC geometricFilter(matchingError, 5000);
  std::vector<size_t> vec_matchingInliers;
  bool valid = geometricFilter.Robust_estimation(featuresI, // points of the query image
                                                 featuresJ, // points of the matched image
                                                 imageSizeI,
                                                 imageSizeJ,
                                                 vec_matchingInliers,
                                                 estimator);
  if(!valid)
  {
    POPART_COUT("[matching]\tGeometric validation failed.");
    return false;
  }
  if(!b_guided_matching)
  {
    // prepare to output vec_featureMatches
    // the indices stored in vec_matchingInliers refer to featuresI, now we need
    // to recover the original indices wrt matchedRegions and queryRegions and fill 
    // a temporary vector.
    std::vector<matching::IndMatch> vec_robustFeatureMatches;
    vec_robustFeatureMatches.reserve(vec_matchingInliers.size());
    for(const size_t idx : vec_matchingInliers)
    {
      // use the index stored in vec_matchingInliers to get the indices of the 
      // original matching features and store them in vec_robustFeatureMatches
      vec_robustFeatureMatches.emplace_back(vec_featureMatches[idx]);
    }
    // just swap the vector so the output is ready
    std::swap(vec_robustFeatureMatches, vec_featureMatches);
  }
  else
  {
    // Use the Fundamental Matrix estimated by the robust estimation to
    // perform guided matching.
    // So we ignore the previous matches and recompute all matches.
    vec_featureMatches.clear();
    geometry_aware::GuidedMatching<
            Mat3,
            openMVG::fundamental::kernel::EpipolarDistanceError>(
            geometricFilter.m_F,
            queryIntrinsicsBase, // cameras::IntrinsicBase of the matched image
            *matcher.getDatabaseRegions(), // features::Regions
            matchedIntrinsicsBase, // cameras::IntrinsicBase of the query image
            matchedRegions._regions, // features::Regions
            Square(geometricFilter.m_dPrecision_robust),
            Square(fDistRatio),
            vec_featureMatches); // output
  }
  return true;
}

bool VoctreeLocalizer::localizeRig(const std::vector<image::Image<unsigned char> > & vec_imageGrey,
                                   const LocalizerParameters *parameters,
                                   std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                                   const std::vector<geometry::Pose3 > &vec_subPoses,
                                   geometry::Pose3 &rigPose, 
                                   std::vector<LocalizationResult> & vec_locResults)
{
  const size_t numCams = vec_imageGrey.size();
  assert(numCams == vec_queryIntrinsics.size());
  assert(numCams == vec_subPoses.size() + 1);

  std::vector<std::unique_ptr<features::Regions> > vec_queryRegions(numCams);
  std::vector<std::pair<std::size_t, std::size_t> > vec_imageSize;
  
  //@todo parallelize?
  for(size_t i = 0; i < numCams; ++i)
  {
    // extract descriptors and features from each image
    vec_queryRegions[i] = std::unique_ptr<features::Regions>(new features::SIFT_Regions());
    POPART_COUT("[features]\tExtract SIFT from query image...");
    _image_describer->Describe(vec_imageGrey[i], vec_queryRegions[i]);
    POPART_COUT("[features]\tExtract SIFT done: found " <<  vec_queryRegions[i]->RegionCount() << " features");
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


bool VoctreeLocalizer::localizeRig(const std::vector<std::unique_ptr<features::Regions> > & vec_queryRegions,
                                   const std::vector<std::pair<std::size_t, std::size_t> > &vec_imageSize,
                                   const LocalizerParameters *parameters,
                                   std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                                   const std::vector<geometry::Pose3 > &vec_subPoses,
                                   geometry::Pose3 &rigPose,
                                   std::vector<LocalizationResult>& vec_locResults)
{
#ifdef HAVE_OPENGV
  if(parameters->_useLocalizeRigNaive)
  {
    POPART_COUT("Using localizeRig_naive()");
    return localizeRig_naive(vec_queryRegions,
                           vec_imageSize,
                           parameters,
                           vec_queryIntrinsics,
                           vec_subPoses,
                           rigPose,
                           vec_locResults);
  }
  else
  {
    POPART_COUT("Using localizeRig_opengv()");
    return localizeRig_opengv(vec_queryRegions,
                              vec_imageSize,
                              parameters,
                              vec_queryIntrinsics,
                              vec_subPoses,
                              rigPose,
                              vec_locResults);
  }
#else
  POPART_COUT("Using localizeRig_naive() -- openGV not supported");
  return localizeRig_naive(vec_queryRegions,
                           vec_imageSize,
                           parameters,
                           vec_queryIntrinsics,
                           vec_subPoses,
                           rigPose,
                           vec_locResults);
#endif
}


#ifdef HAVE_OPENGV

bool VoctreeLocalizer::localizeRig_opengv(const std::vector<std::unique_ptr<features::Regions> > & vec_queryRegions,
                                          const std::vector<std::pair<std::size_t, std::size_t> > &vec_imageSize,
                                          const LocalizerParameters *parameters,
                                          std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                                          const std::vector<geometry::Pose3 > &vec_subPoses,
                                          geometry::Pose3 &rigPose,
                                          std::vector<LocalizationResult>& vec_locResults)
{
  const std::size_t numCams = vec_queryRegions.size();
  assert(numCams == vec_queryIntrinsics.size());
  assert(numCams == vec_subPoses.size() + 1);   
  
  const VoctreeLocalizer::Parameters *param = static_cast<const VoctreeLocalizer::Parameters *>(parameters);
  if(!param)
  {
    // error!
    throw std::invalid_argument("The parameters are not in the right format!!");
  }

  std::vector<std::map< pair<IndexT, IndexT>, std::size_t > > vec_occurrences(numCams);
  std::vector<std::vector<voctree::DocMatch> > vec_matchedImages(numCams);
  std::vector<Mat> vec_pts3D(numCams);
  std::vector<Mat> vec_pts2D(numCams);

  // for each camera retrieve the associations
  //@todo parallelize?
  size_t numAssociations = 0;
  for(std::size_t cam = 0; cam < numCams; ++cam)
  {

    // this map is used to collect the 2d-3d associations as we go through the images
    // the key is a pair <Id3D, Id2d>
    // the element is the pair 3D point - 2D point
    auto &occurrences = vec_occurrences[cam];
    auto &matchedImages = vec_matchedImages[cam];
    auto &imageSize = vec_imageSize[cam];
    Mat &pts3D = vec_pts3D[cam];
    Mat &pts2D = vec_pts2D[cam];
    cameras::Pinhole_Intrinsic_Radial_K3 &queryIntrinsics = vec_queryIntrinsics[cam];
    features::SIFT_Regions &queryRegions = *dynamic_cast<features::SIFT_Regions*> (vec_queryRegions[cam].get());
    const bool useInputIntrinsics = true;
    getAllAssociations(queryRegions,
                       imageSize,
                       *param,
                       useInputIntrinsics,
                       queryIntrinsics,
                       occurrences,
                       pts2D,
                       pts3D,
                       matchedImages);
    numAssociations += occurrences.size();
  }
  
  // @todo Here it could be possible to filter the associations according to their
  // occurrences, eg giving priority to those associations that are more frequent

  const std::size_t minNumAssociations = 5;  //possible parameter?
  if(numAssociations < minNumAssociations)
  {
    POPART_COUT("[poseEstimation]\tonly " << numAssociations << " have been found, not enough to do the resection!");
    return false;
  }
  
  {
    for(std::size_t cam = 0; cam < vec_subPoses.size(); ++cam)
    {
    POPART_COUT("Rotation: " << vec_subPoses[cam].rotation());
    POPART_COUT("Centre: " << vec_subPoses[cam].center());
    }
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
  
  { // just debugging stuff, this block can be removed
    
    if(vec_inliers.size() < numCams)
    {
      // in general the inlier should be spread among different cameras
      POPART_CERR("WARNING: RIG Voctree Localizer: Inliers in " 
              << vec_inliers.size() << " cameras on a RIG of " 
              << numCams << " cameras.");
    }

    for(std::size_t cam = 0; cam < vec_inliers.size(); ++cam)
      POPART_COUT("#inliers for cam " << cam << ": " << vec_inliers[cam].size());
    
    POPART_COUT("Pose after resection:");
    POPART_COUT("Rotation: " << rigPose.rotation());
    POPART_COUT("Centre: " << rigPose.center());
    
    // compute the reprojection error for inliers (just debugging purposes)
    for(std::size_t cam = 0; cam < numCams; ++cam)
    {
      const std::size_t numPts = vec_pts2D[cam].cols();
      const cameras::Pinhole_Intrinsic_Radial_K3 &currCamera = vec_queryIntrinsics[cam];
      Mat2X residuals;
      if(cam!=0)
        residuals = currCamera.residuals(vec_subPoses[cam-1]*rigPose, vec_pts3D[cam], vec_pts2D[cam]);
      else
        residuals = currCamera.residuals(geometry::Pose3()*rigPose, vec_pts3D[cam], vec_pts2D[cam]);

      auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();

//      POPART_COUT("Camera " << cam << " all reprojection errors:");
//      POPART_COUT(sqrErrors);
//
//      POPART_COUT("Camera " << cam << " inliers reprojection errors:");
      const auto &currInliers = vec_inliers[cam];

      double rmse = 0;
      for(std::size_t j = 0; j < currInliers.size(); ++j)
      {
//          std::cout << sqrErrors(currInliers[j]) << " ";
          rmse += sqrErrors(currInliers[j]);
      }
      if(!currInliers.empty())
        POPART_COUT("\n\nRMSE inliers: " << std::sqrt(rmse/currInliers.size()));
    }
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
      pose = vec_subPoses[cam-1] * rigPose;
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
    
    const auto &matchedImages = vec_matchedImages[cam];
    
    vec_locResults.emplace_back(matchData, indMatch3D2D, pose, intrinsics, matchedImages, refineOk);
  }
  
  
  if(!refineOk)
    return resectionOk;
  
  { // just debugging stuff, this block can be removed
    
    POPART_COUT("Pose after BA:");
    POPART_COUT("Rotation: " << rigPose.rotation());
    POPART_COUT("Centre: " << rigPose.center());
    
    // compute the reprojection error for inliers (just debugging purposes)
    for(std::size_t cam = 0; cam < numCams; ++cam)
    {
      const std::size_t numPts = vec_pts2D[cam].cols();
      const cameras::Pinhole_Intrinsic_Radial_K3 &currCamera = vec_queryIntrinsics[cam];
      Mat2X residuals;
      if(cam!=0)
        residuals = currCamera.residuals(vec_subPoses[cam-1]*rigPose, vec_pts3D[cam], vec_pts2D[cam]);
      else
        residuals = currCamera.residuals(geometry::Pose3()*rigPose, vec_pts3D[cam], vec_pts2D[cam]);

      auto sqrErrors = (residuals.cwiseProduct(residuals)).colwise().sum();

//      POPART_COUT("Camera " << cam << " all reprojection errors after BA:");
//      POPART_COUT(sqrErrors);

//      POPART_COUT("Camera " << cam << " inliers reprojection errors after BA:");
      const auto &currInliers = vec_inliers[cam];

      double rmse = 0;
      for(std::size_t j = 0; j < currInliers.size(); ++j)
      {
//          std::cout << sqrErrors(currInliers[j]) << " ";
          rmse += sqrErrors(currInliers[j]);
      }
      if(!currInliers.empty())
        POPART_COUT("\n\nRMSE inliers cam " << cam << ": " << std::sqrt(rmse/currInliers.size()));
    }
  }
    
  return resectionOk;
}

#endif // HAVE_OPENGV

// subposes is n-1 as we consider the first camera as the main camera and the 
// reference frame of the grid
bool VoctreeLocalizer::localizeRig_naive(const std::vector<std::unique_ptr<features::Regions> > & vec_queryRegions,
                                          const std::vector<std::pair<std::size_t, std::size_t> > &vec_imageSize,
                                          const LocalizerParameters *parameters,
                                          std::vector<cameras::Pinhole_Intrinsic_Radial_K3 > &vec_queryIntrinsics,
                                          const std::vector<geometry::Pose3 > &vec_subPoses,
                                          geometry::Pose3 &rigPose,
                                          std::vector<LocalizationResult>& vec_localizationResults)
{
  const size_t numCams = vec_queryRegions.size();
  
  assert(numCams==vec_queryIntrinsics.size());
  assert(numCams==vec_subPoses.size()+1);
  assert(numCams==vec_imageSize.size());

  vec_localizationResults.resize(numCams);
    
  // this is basic, just localize each camera alone
  //@todo parallelize?
  std::vector<bool> isLocalized(numCams, false);
  for(size_t i = 0; i < numCams; ++i)
  {
    isLocalized[i] = localize(vec_queryRegions[i], vec_imageSize[i], parameters, true /*useInputIntrinsics*/, vec_queryIntrinsics[i], vec_localizationResults[i]);
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
    return vec_localizationResults[0].isValid();
  }
  
  // if only one camera has been localized
  if(numLocalizedCam == 1)
  {
    // all the other cameras have not been localized just return the result of the 
    // localized one
    
    // find the index of the localized camera
    const std::size_t idx = std::distance(isLocalized.begin(), 
                                          std::find(isLocalized.begin(), isLocalized.end(), true));
    
    // useless safeguard as there should be at least 1 element at this point but
    // better safe than sorry
    assert(idx < isLocalized.size());
    
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
    return vec_localizationResults[idx].isValid();
  }
  
  // ** otherwise run a BA with the localized cameras

  refineRigPose(vec_subPoses, vec_localizationResults, rigPose);
  
  return true;
}



} // localization
} // openMVG
