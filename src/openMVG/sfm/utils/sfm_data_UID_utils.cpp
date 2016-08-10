#include "sfm_data_UID_utils.hpp"
#include "openMVG/exif/exif_IO_EasyExif.hpp"

namespace openMVG {
namespace sfm {

void updateStructureWithNewUID(Landmarks &landmarks, const std::map<std::size_t, std::size_t> &oldIdToNew)
{
  // update the id in the visibility of each 3D point
  for(auto &iter : landmarks)
  {
    Landmark& currentLandmark = iter.second;
    
    // the new observations where to copy the existing ones
    // (needed as the key of the map is the idview)
    Observations newObservations;
    
    for(const auto &iterObs : currentLandmark.obs)
    {
      const auto idview = iterObs.first;
      const Observation &obs = iterObs.second;

      newObservations.emplace(oldIdToNew.at(idview), obs);
    }
    
    assert(currentLandmark.obs.size() == newObservations.size());
    currentLandmark.obs.swap(newObservations);
  }  
}


void sanityCheckLandmarks(const Landmarks &landmarks, const Views &views)
{
  for(const auto &iter : landmarks)
  {
    const Landmark& currentLandmark = iter.second;
    for(const auto &iterObs : currentLandmark.obs)
    {
      const auto idview = iterObs.first;
      const Observation &obs = iterObs.second;

      // there must be a view with that id (in the map) and the view must have 
      // the same id (the member)
      assert(views.count(idview) == 1);
      assert(views.at(idview)->id_view == idview);
    }
  }  
}

void regenerateUID(SfM_Data &sfmdata, std::map<std::size_t, std::size_t> &oldIdToNew, bool sanityCheck)
{
  // if the views are empty, nothing to be done. 
  if(sfmdata.GetViews().empty())
    return;
  
  updateViewIDs(sfmdata.views, oldIdToNew);
  
  if(!sanityCheck)
    return;
  
  sanityCheckLandmarks(sfmdata.GetLandmarks(), sfmdata.GetViews());
  
  sanityCheckLandmarks(sfmdata.GetControl_Points(), sfmdata.GetViews());
  
}


void updateViewIDs(Views &views, std::map<std::size_t, std::size_t> &oldIdToNew)
{
  // if the views are empty, nothing to be done. 
  if(views.empty())
    return;
  
  Views newViews;

  for(auto const &iter : views)
  {
    const View& currentView = *iter.second.get();
    const auto &imageName = currentView.s_Img_path;
    
    exif::Exif_IO_EasyExif exifReader(imageName);

    // compute the view UID
    const std::size_t uid = exif::computeUID(exifReader, imageName);

    // update the mapping
    assert(oldIdToNew.count(currentView.id_view) == 0);
    oldIdToNew.emplace(currentView.id_view, uid);
    
    // add the view to the new map using the uid as key and change the id
    assert(newViews.count(uid) == 0);
    newViews.emplace(uid, iter.second);
    newViews[uid]->id_view = uid;
  }
  
  assert(newViews.size() == views.size());
  views.swap(newViews);
}

}
}