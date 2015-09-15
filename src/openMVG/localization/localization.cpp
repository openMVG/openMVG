#include "localization.hpp"


#include <openMVG/sfm/sfm_data_io.hpp>
#include <openMVG/sfm/pipelines/sfm_robust_model_estimation.hpp>
#include <openMVG/features/io_regions_type.hpp>
#include <openMVG/matching/regions_matcher.hpp>
#include <openMVG/matching_image_collection/Matcher.hpp>
#include <openMVG/matching/matcher_kdtree_flann.hpp>
#include <openMVG/robust_estimation/guided_matching.hpp>
#include <openMVG/matching_image_collection/F_ACRobust.hpp>
#include <third_party/progress/progress.hpp>
//#include <cereal/archives/json.hpp>

#include <algorithm>

//@fixme move/redefine
#define POPART_COUT(x) std::cout << x << std::endl
#define POPART_CERR(x) std::cerr << x << std::endl

namespace openMVG {
namespace localization {


//// maybe better as a static method of Image_describer
//bool Init_image_describer_type_from_file(const std::string & sImage_describer, 
//                                         features::Image_describer & image_describer)
//{
//  if (!stlplus::is_file(sImage_describer))
//  {
//    std::cerr << "Expected file image_describer.json cannot be opened." << std::endl;
//    return false;
//  }
//
//  // Dynamically load the image_describer from the file (will restore old used settings)
//  std::ifstream stream(sImage_describer.c_str());
//  if (!stream.is_open())
//  {
//    std::cerr << "Unable to open image_describer.json." << std::endl;
//    return false;
//  }
//
//  try
//  {
//    cereal::JSONInputArchive archive(stream);
//    archive(cereal::make_nvp("image_describer", image_describer));
//    return true;
//  }
//  catch (const cereal::Exception & e)
//  {
//    std::cerr << e.what() << std::endl
//      << "Cannot dynamically allocate the Image_describer interface." << std::endl;
//    return false;
//  }
//}



// inputs
// - sfmdata path
// - descriptorsFolder directory with the sift
// - vocTreeFilepath; 
// - weightsFilepath; 
bool VoctreeLocalizer::init( const std::string &sfmFilePath,
                            const std::string &descriptorsFolder,
                            const std::string &vocTreeFilepath,
                            const std::string &weightsFilepath)
{
  using namespace openMVG::features;
  
  // load the sfm data containing the 3D reconstruction info
  if (!Load(_sfm_data, sfmFilePath, sfm::ESfM_Data::ALL)) 
  {
    std::cerr << std::endl
      << "The input SfM_Data file "<< sfmFilePath << " cannot be read." << std::endl;
    return false;
  }

  // load the features and descriptors
  // initially we need all the feature in order to create the database
  // then we can store only those associated to 3D points
  //? can we use Feature_Provider to load the features and filter them later?

  // this block is used to get the type of features (by default SIFT) used
  // for the reconstruction
  const std::string sImage_describer = stlplus::create_filespec(descriptorsFolder, "image_describer", "json");
  std::unique_ptr<Regions> regions_type = Init_region_type_from_file(sImage_describer);
  if(!regions_type)
  {
    std::cerr << "Invalid: "
            << sImage_describer << " regions type file." << std::endl;
    return false;
  }

//  // initialize the image describer from image_describer.json in matches directory
//  // @fixme it would be maybe better to initialiaze to SIFT by default if it fails
//  // to garantee the back compatibility
//  if(!Init_image_describer_type_from_file(sImage_describer, *_image_describer.get()))
//  {
//    std::cerr << "Unable to initialize the image describer." << std::endl;
//    return false;
//  }
  
//  // Load all the features and store them inside the region provider
//  std::shared_ptr<sfm::Regions_Provider> regions_provider = std::make_shared<sfm::Regions_Provider>();
//  if(!regions_provider->load(_sfm_data, descriptorsFolder, regions_type)) // instead of passing by Init_region_type_from_file to get the type
//  {
//    std::cerr << std::endl << "Invalid regions." << std::endl;
//    return false;
//  }
    
  initDatabase(vocTreeFilepath, weightsFilepath, descriptorsFolder);

  // filter the features to keep only those having a 3D point
  // associated to them. loadReconstructionDescriptors() but just with the 
  // filtering part
  
  return true;
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
  C_Progress_display my_progress_bar(_sfm_data.GetViews().size(),
                                     std::cout, "\n- Regions Loading -\n");

  POPART_COUT("Build observations per view");
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

  POPART_COUT("Load Features and Descriptors per view");
  // Read for each view the corresponding regions and store them
  for(const auto &iter : _sfm_data.GetViews())
  {
    const std::shared_ptr<sfm::View> currView = iter.second;
    const IndexT id_view = currView->id_view;
    Reconstructed_RegionsT& currRecoRegions = _regions_per_view[id_view];

    const std::string sImageName = stlplus::create_filespec(_sfm_data.s_root_path, currView.get()->s_Img_path);
    const std::string basename = stlplus::basename_part(sImageName);
    const std::string featFilepath = stlplus::create_filespec(feat_directory, basename, ".feat");
    const std::string descFilepath = stlplus::create_filespec(feat_directory, basename, ".desc");
    //    std::cout << "Feat: " << featFilepath << std::endl;
    //    std::cout << "Desc: " << descFilepath << std::endl;

    if(!currRecoRegions._regions.Load(featFilepath, descFilepath))
    {
      std::cerr << "Invalid regions files for the view: " << sImageName << std::endl;
      return false;
    }
    
    std::vector<voctree::Word> words = _voctree.quantize(currRecoRegions._regions.Descriptors());
    voctree::DocId docId = _database.insert(words);
    _mapDocIdToView[docId] = id_view;

    // Filter descriptors to keep only the 3D reconstructed points
    currRecoRegions.filterRegions(observationsPerView[id_view]);
    ++my_progress_bar;
  }
  
  return true;
}


bool VoctreeLocalizer::Localize( const image::Image<unsigned char> & imageGray,
                cameras::IntrinsicBase * queryIntrinsics,
                const size_t numResults,
                geometry::Pose3 & pose,
                bool useGuidedMatching,
                sfm::Image_Localizer_Match_Data * resection_data /*= nullptr*/)
{
  // extract descriptors and features from image
  POPART_COUT("[features]\tExtract SIFT from query image");
  std::unique_ptr<features::Regions> tmpQueryRegions(new features::SIFT_Regions());
  _image_describer.Describe(imageGray, tmpQueryRegions, NULL);
  POPART_COUT("[features]\tExtract SIFT done: found " << tmpQueryRegions->RegionCount() << " features");
  features::SIFT_Regions queryRegions = *dynamic_cast<features::SIFT_Regions*> (tmpQueryRegions.get());
//  POPART_COUT("Extract SIFT done!!");

  POPART_COUT("[database]\tRequest closest images from voctree");
  // pass the descriptors through the vocabulary tree to get the visual words
  // associated to each feature
  std::vector<voctree::Word> requestImageWords = _voctree.quantize(queryRegions.Descriptors());
  
  // Request closest images from voctree
  std::vector<voctree::Match> matchedImages;
  _database.find(requestImageWords, numResults, matchedImages);
  
  // just debugging bla bla
  // for each similar image found print score and number of features
  for(const voctree::Match & currMatch : matchedImages )
  {
    // get the corresponding index of the view
    const IndexT matchedViewIndex = _mapDocIdToView[currMatch.id];
    // get the view handle
    const std::shared_ptr<sfm::View> matchedView = _sfm_data.views[matchedViewIndex];
    POPART_COUT( "[database]\t\t match " << matchedView->s_Img_path 
            << " [docid: "<< currMatch.id << "]"
            << " with score " << currMatch.score 
            << " and it has "  << _regions_per_view[matchedViewIndex]._regions.RegionCount() 
            << " features with 3D points");
  }

  //@fixme Maybe useless, just do everything with DistanceRatioMatch
  // preparing the matcher, it will use the extracted Regions as reference and it
  // will match them to the Regions of each similar image
  const float fDistRatio = 0.6; //@fixme this could be a param
  typedef flann::L2<unsigned char> MetricT;
  typedef matching::ArrayMatcher_Kdtree_Flann<unsigned char, MetricT> MatcherT;
  POPART_COUT("[matching]\tBuild the matcher");
  matching::RegionsMatcherT<MatcherT> matcher(queryRegions);
  
  // Prepare intrinsics @fixme 
  POPART_COUT("[matching]\tPrepare query intrinsics");
  const bool bKnownIntrinsic = (queryIntrinsics != nullptr);
  Mat3 K = Mat3::Identity();
 
  // for each found similar image
  for(const voctree::Match& matchedImage : matchedImages)
  {
    // the view index of the current matched image
    const IndexT matchedViewIndex = _mapDocIdToView[matchedImage.id];
    // the handler to the current view
    const std::shared_ptr<sfm::View> matchedView = _sfm_data.views[matchedViewIndex];
    // its associated reconstructed regions
    const Reconstructed_RegionsT& matchedRegions = _regions_per_view[matchedViewIndex];
    // its associated intrinsics
    const cameras::IntrinsicBase *matchedIntrinsics = _sfm_data.intrinsics[matchedView->id_intrinsic].get();
    std::vector<matching::IndMatch> vec_featureMatches;
    bool matchWorked = robustMatching( matcher, 
                                      queryIntrinsics,
                                      matchedRegions,
                                      matchedIntrinsics,
                                      fDistRatio,
                                      useGuidedMatching,
                                      std::make_pair(imageGray.Width(), imageGray.Height()),
                                      std::make_pair(imageGray.Width(), imageGray.Height()), // NO! @fixme here we need the size of the img in the dataset...
                                      vec_featureMatches);
    if (!matchWorked)
    {
      POPART_COUT("\tmatching with " << matchedView->s_Img_path << " failed! Skipping image");
      continue;
    }
  
    // recover the 2D-3D associations
    // Prepare data for resection
    sfm::Image_Localizer_Match_Data matchData;
    matchData.pt2D = Mat2X(2, vec_featureMatches.size());
    matchData.pt3D = Mat3X(3, vec_featureMatches.size());

    // Get the 3D points associated to each matched feature
    std::size_t index = 0;
    for(const matching::IndMatch& featureMatch : vec_featureMatches)
    {
      // the ID of the 3D point
      IndexT trackId3D = matchedRegions._associated3dPoint[featureMatch._j];

      // prepare data for resectioning
      matchData.pt3D.col(index) = _sfm_data.GetLandmarks().at(trackId3D).X;

      const Vec2 feat = queryRegions.GetRegionPosition(featureMatch._i);
      if(bKnownIntrinsic)
        matchData.pt2D.col(index) = queryIntrinsics->get_ud_pixel(feat);
      else
        matchData.pt2D.col(index) = feat;

      ++index;
    }
    // estimate the pose
    // Do the resectioning: compute the camera pose.
    std::vector<size_t> vec_inliers;
    double errorMax = std::numeric_limits<double>::max();

    bool bResection = sfm::robustResection(
                                           std::make_pair(imageGray.Width(), imageGray.Height()),
                                           matchData.pt2D, matchData.pt3D,
                                           &vec_inliers,
                                           // Use intrinsic guess if possible
                                           (bKnownIntrinsic) ? &K : NULL,
                                           &matchData.projection_matrix, &errorMax);

    std::cout << std::endl
            << "-------------------------------" << std::endl
            << "-- Robust Resection using view: " << _mapDocIdToView[matchedImage.id] << std::endl
            << "-- Resection status: " << bResection << std::endl
            << "-- #Points used for Resection: " << vec_featureMatches.size() << std::endl
            << "-- #Points validated by robust Resection: " << vec_inliers.size() << std::endl
            << "-- Threshold: " << errorMax << std::endl
            << "-------------------------------" << std::endl;

    if(!bResection)
    {
      POPART_COUT("\tResection FAILED");
      continue;
    }
    POPART_COUT("Resection SUCCEDED");
    // Decompose P matrix
    Mat3 K_, R_;
    Vec3 t_;
    KRt_From_P(matchData.projection_matrix, &K_, &R_, &t_);
    pose = geometry::Pose3(R_, -R_.transpose() * t_);
    
    bool refineStatus = sfm::SfM_Localizer::RefinePose(queryIntrinsics, pose, matchData, true /*b_refine_pose*/, false /*b_refine_intrinsic*/);
    if(!refineStatus)
      POPART_COUT("Refine pose could not improve the estimation of the camera pose.");
    POPART_COUT("K: " << K_);
    break;
  }

  return true;
  
 } 
  





bool VoctreeLocalizer::robustMatching(matching::RegionsMatcherT<MatcherT> & matcher, 
                    const cameras::IntrinsicBase * queryIntrinsics,   // the intrinsics of the image we are using as reference
                    const Reconstructed_RegionsT & matchedRegions,
                    const cameras::IntrinsicBase * matchedIntrinsics,
                    const float fDistRatio,
                    const bool b_guided_matching,
                    const std::pair<size_t,size_t> & imageSizeI,     // size of the first image @fixme change the API of the kernel!! 
                    const std::pair<size_t,size_t> & imageSizeJ,     // size of the first image
                    std::vector<matching::IndMatch> & vec_featureMatches) const
{
  // A. Putative Features Matching
  bool matchWorked = matcher.Match(fDistRatio, matchedRegions._regions, vec_featureMatches);
  if (!matchWorked)
  {
    POPART_COUT("\tmatching with failed!");
    return false;
  }

  // they contains the matching features
  Mat featuresI(2, vec_featureMatches.size());
  Mat featuresJ(2, vec_featureMatches.size());

  // fill the matrices with the features according to vec_featureMatches
  for(int i = 0; i < vec_featureMatches.size(); ++i)
  {
    const matching::IndMatch& match = vec_featureMatches[i];
    //@fixme get the query regions, we need to add the method to RegionsMatcherT
    featuresI.col(i) = matcher.getDatabaseRegions()->GetRegionPosition(match._i);
    featuresJ.col(i) = matchedRegions._regions.GetRegionPosition(match._j);
  }

  matching_image_collection::GeometricFilter_FMatrix_AC geometricFilter(4.0);
  std::vector<size_t> vec_matchingInliers;
  bool valid = geometricFilter.Robust_estimation(featuresI, // points of the query image
                                                 featuresJ, // points of the matched image
                                                 imageSizeI,
                                                 imageSizeJ,
                                                 vec_matchingInliers);
  if(!valid)
  {
    POPART_COUT("\tUnable to robustly matching the query image with the database image.");
    return false;
  }
  if(!b_guided_matching)
  {
    std::vector<matching::IndMatch> vec_robustFeatureMatches(vec_matchingInliers.size());
    for(const int i : vec_matchingInliers)
    {
      vec_robustFeatureMatches[i] = vec_featureMatches[i];
    }
    // replace the featuresMatches with the robust ones.
    std::swap(vec_featureMatches, vec_robustFeatureMatches);
  }
  else
  {
    // Use the Fundamental Matrix estimated by the robust estimation to
    // perform guided matching.
    // So we ignore the previous matches and recompute all matches.

    geometry_aware::GuidedMatching<
            Mat3,
            openMVG::fundamental::kernel::EpipolarDistanceError>(
            geometricFilter.m_F,
            queryIntrinsics, // cameras::IntrinsicBase of the matched image
            *matcher.getDatabaseRegions(), // features::Regions
            matchedIntrinsics, // cameras::IntrinsicBase of the query image
            matchedRegions._regions, // features::Regions
            Square(geometricFilter.m_dPrecision_robust),
            Square(fDistRatio),
            vec_featureMatches); // output
  }
  
  return true;
}

} // localization
} // openMVG
