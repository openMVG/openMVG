#include "localization.hpp"

#include "third_party/progress/progress.hpp"

#include <algorithm>


namespace openMVG {
namespace localization {


bool Localization::loadReconstructionDescriptors(
  const sfm::SfM_Data & sfm_data,
  const std::string & feat_directory)
{
  C_Progress_display my_progress_bar( sfm_data.GetViews().size(),
    std::cout, "\n- Regions Loading -\n");
  
  std::cout << "Build observations per view" << std::endl;
  // Build observations per view
  std::map<IndexT, std::vector<FeatureInImage> > observationsPerView;
  for( auto landmarkValue: sfm_data.structure )
  {
    IndexT trackId = landmarkValue.first;
    sfm::Landmark& landmark = landmarkValue.second;
    for(auto obs: landmark.obs)
    {
      const IndexT viewId = obs.first;
      const sfm::Observation& obs2d = obs.second;
      observationsPerView[viewId].push_back(FeatureInImage(obs2d.id_feat, trackId));
    }
  }
  for( auto featuresInImage: observationsPerView )
  {
    std::sort(featuresInImage.second.begin(), featuresInImage.second.end());
  }

  std::cout << "Load Features and Descriptors per view" << std::endl;
  // Read for each view the corresponding regions and store them
  for (sfm::Views::const_iterator iter = sfm_data.GetViews().begin();
    iter != sfm_data.GetViews().end(); ++iter, ++my_progress_bar)
  {
    const IndexT id_view = iter->second->id_view;
    Reconstructed_RegionsT& reconstructedRegion = regions_per_view[id_view];

    const std::string sImageName = stlplus::create_filespec(sfm_data.s_root_path, iter->second.get()->s_Img_path);
    const std::string basename = stlplus::basename_part(sImageName);
    const std::string featFilepath = stlplus::create_filespec(feat_directory, basename, ".feat");
    const std::string descFilepath = stlplus::create_filespec(feat_directory, basename, ".desc");
//    std::cout << "Feat: " << featFilepath << std::endl;
//    std::cout << "Desc: " << descFilepath << std::endl;

    if (!reconstructedRegion._regions.Load(featFilepath, descFilepath))
    {
      std::cerr << "Invalid regions files for the view: " << sImageName << std::endl;
      return false;
    }

    // Filter descriptors to only keep the 3D reconstructed points
    reconstructedRegion.filterRegions(observationsPerView[id_view]);
  }
  return true;
}

}
}
