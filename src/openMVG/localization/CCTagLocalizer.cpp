#include "CCTagLocalizer.hpp"
#include "openMVG/sfm/sfm_data.hpp"
#include "reconstructed_regions.hpp"
#include <openMVG/sfm/sfm_data_io.hpp>

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
  
  for(const auto & landmark : landmarks)
  {
    // Use the first observation to retrieve the associated descriptor.
    const auto & firstObservation = *landmark.second.obs.begin();
    
    // Retrieve the Regions of the first observation
    auto & reconstructedRegions = _regions_per_view[firstObservation.first];
    
    // Get the feature id: remap the index as we only load the reconstructed regions
    const auto localFeatureId = reconstructedRegions._mapFullToLocal[firstObservation.second.id_feat];
    
    const auto & desc = reconstructedRegions._regions.Descriptors()[localFeatureId];
    IndexT idCCTag = getCCTagId(desc);

    // Insert <idCCTag, 3D point> into a map.
    if (idCCTag!=UndefinedIndexT)
    {
      _cctagDatabase.emplace(idCCTag, landmark.second.X);
    }
  }
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
                cameras::Pinhole_Intrinsic &queryIntrinsics,
                geometry::Pose3 & pose,
                sfm::Image_Localizer_Match_Data &resection_data,
                std::vector<pair<IndexT, IndexT> > &associationIDs)
{
#if 0
  // extract descriptors and features from image
  POPART_COUT("[features]\tExtract CCTag from query image");
  std::unique_ptr<features::Regions> tmpQueryRegions(new features::CCTAG_Regions());
  _image_describer.Describe(imageGray, tmpQueryRegions, NULL);
  POPART_COUT("[features]\tExtract CCTAG done: found " << tmpQueryRegions->RegionCount() << " features");
  features::CCTAG_Regions queryRegions = *dynamic_cast<features::CCTAG_Regions*> (tmpQueryRegions.get());
//  
//    // recover the 2D-3D associations
//    // Prepare data for resection
  
  // Mapping between cctag id and their index inside Regions
  std::map<IndexT, IndexT> cctagWith3D;
  
  for(size_t i = 0; i<queryRegions.RegionCount(); ++i)
  {
    // get the current descriptor
    const CCTagDescriptor & desc = queryRegions.Descriptors()[i];
    IndexT idCCTag = getCCTagId(desc);
    
    // check whether it is in the database of the cctag with an associated 3D point
    if(_cctagDatabase.find(idCCTag) != _cctagDatabase.end())
    {
      cctagWith3D[idCCTag] = i;
    }
  }
  
    sfm::Image_Localizer_Match_Data matchData;
    matchData.pt2D = Mat2X(2, cctagWith3D.size());
    matchData.pt3D = Mat3X(3, cctagWith3D.size());
#endif
    
//
//    // Get the 3D points associated to each matched feature
//    std::size_t index = 0;
//    for(const matching::IndMatch& featureMatch : vec_featureMatches)
//    {
//      // the ID of the 3D point
//      IndexT trackId3D = matchedRegions._associated3dPoint[featureMatch._j];
//
//      // prepare data for resectioning
//      matchData.pt3D.col(index) = _sfm_data.GetLandmarks().at(trackId3D).X;
//
//      const Vec2 feat = queryRegions.GetRegionPosition(featureMatch._i);
//      if(bKnownIntrinsic)
//        matchData.pt2D.col(index) = queryIntrinsics->get_ud_pixel(feat);
//      else
//        matchData.pt2D.col(index) = feat;
//
//      ++index;
//    }
//    // estimate the pose
//    // Do the resectioning: compute the camera pose.
//    std::vector<size_t> vec_inliers;
//    double errorMax = std::numeric_limits<double>::max();
//
//    bool bResection = sfm::robustResection(std::make_pair(imageGray.Width(), imageGray.Height()),
//                                           matchData.pt2D, matchData.pt3D,
//                                           &vec_inliers,
//                                           // Use intrinsic guess if possible
//                                           (bKnownIntrinsic) ? &K : NULL,
//                                           &matchData.projection_matrix, &errorMax);
//
//    std::cout << std::endl
//            << "-------------------------------" << std::endl
//            << "-- Robust Resection using view: " << _mapDocIdToView[matchedImage.id] << std::endl
//            << "-- Resection status: " << bResection << std::endl
//            << "-- #Points used for Resection: " << vec_featureMatches.size() << std::endl
//            << "-- #Points validated by robust Resection: " << vec_inliers.size() << std::endl
//            << "-- Threshold: " << errorMax << std::endl
//            << "-------------------------------" << std::endl;
//
//    if(!bResection)
//    {
//      POPART_COUT("\tResection FAILED");
//      continue;
//    }
//    POPART_COUT("Resection SUCCEDED");
//    // Decompose P matrix
//    Mat3 K_, R_;
//    Vec3 t_;
//    KRt_From_P(matchData.projection_matrix, &K_, &R_, &t_);
//    pose = geometry::Pose3(R_, -R_.transpose() * t_);
//    
//    bool refineStatus = sfm::SfM_Localizer::RefinePose(queryIntrinsics, pose, matchData, true /*b_refine_pose*/, false /*b_refine_intrinsic*/);
//    if(!refineStatus)
//      POPART_COUT("Refine pose could not improve the estimation of the camera pose.");
//    POPART_COUT("K: " << K_);
//    break;
//  }
//
return true;
  
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

} // localization
} // openMVG
