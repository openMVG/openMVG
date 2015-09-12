#include "localization.hpp"


#include <openMVG/sfm/sfm_data_io.hpp>
#include <openMVG/features/io_regions_type.hpp>
#include <third_party/progress/progress.hpp>
#include <cereal/archives/json.hpp>

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
    Reconstructed_RegionsT& reconstructedRegion = _regions_per_view[id_view];

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


}
}
